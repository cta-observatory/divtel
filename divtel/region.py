"""
A weighted region of sky, and what an array sees of it.

A pointing strategy needs a target, and for anything but a point source the
target is a *region* with structure: a gravitational-wave localization is a
probability distribution over the sky, a survey field is a patch with edges, a
galaxy catalogue is a list of directions with weights on them. `SkyRegion` is
that target reduced to the only two things the geometry needs -- a set of unit
vectors and a weight on each -- and everything in `divtel.strategy` works on it.

Reducing it that far is what keeps this module free of any sky-map format.
Building a region from a published HEALPix localization needs a HEALPix reader
and lives in `divtel.skymap`, behind an optional dependency; once built, nothing
downstream knows or cares where the directions came from.

Frames
------
The vectors are unit vectors in whatever frame they were built in, and two
frames are used here for two different questions.

`from_coords` transforms sky coordinates into the array's horizontal frame at a
given time, which answers "what could the array have done that night". The
region sits where the sky actually put it.

`from_table` and `from_lonlat` keep the region in its own celestial frame, and
`place` then rotates it rigidly to a chosen altitude and azimuth, which answers
"what can an array do with a region of this shape". The elevation has to be
controllable for that question, because an array looking near the horizon is
foreshortened and the same divergence buys a different spread at 20 degrees
than at 70.

Sizes
-----
"How big is it" has three answers and the gap between them is the point. A
two-detector localization is a long thin arc, so its area and its length say
completely different things: `equivalent_radius` treats it as a disc,
`containment_radius` gives the cone an array's fields of view have to span, and
`max_separation` gives the region's length. For GW170817's first alert map those
are 7.7, 62.1 and 124.3 degrees.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field, replace
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

from .pointing import alt_az_to_vector, as_altaz

__all__ = ["SkyRegion", "as_altaz", "rotation_between"]

# Above this many directions the exact pairwise sweep in `max_separation` is
# replaced by an evenly-spaced subsample. At the pixel scale of a real sky map
# the two agree to far better than the accuracy anyone needs from a diameter.
_MAX_EXACT_PIXELS = 20_000


def rotation_between(source, target):
    """
    The shortest rotation carrying one direction onto another.

    Parameters
    ----------
    source, target: `numpy.ndarray`
        unit vectors

    Returns
    -------
    `numpy.ndarray`
        a 3 by 3 rotation matrix, to apply as ``vectors @ rotation.T``
    """
    axis = np.cross(source, target)
    sine = float(np.linalg.norm(axis))
    cosine = float(np.dot(source, target))

    if sine < 1e-12:
        return np.eye(3) if cosine > 0 else -np.eye(3)

    axis = axis / sine
    cross = np.array([[0.0, -axis[2], axis[1]],
                      [axis[2], 0.0, -axis[0]],
                      [-axis[1], axis[0], 0.0]])
    angle = np.arctan2(sine, cosine)
    return (np.eye(3) + np.sin(angle) * cross
            + (1 - np.cos(angle)) * (cross @ cross))


@dataclass(frozen=True)
class SkyRegion:
    """
    A set of sky directions with a weight on each.

    Attributes
    ----------
    directions: `numpy.ndarray`
        shape (n, 3), unit vectors. See the module docstring on which frame.
    weights: `numpy.ndarray`
        shape (n,), the weight of each direction -- a probability, for a
        localization. Not normalised: `visible_part` cuts directions away and
        leaves the survivors' weights alone, so that coverage measured
        afterwards stays a fraction of the whole original region rather than
        silently renormalising the abandoned part out of existence.
    pixel_area: `astropy.Quantity` or None
        solid angle each direction stands for, if the region came from a
        pixelization. Without it `area` and `equivalent_radius` are unavailable.
    meta: dict
        provenance -- which map, which credible level, which circular carried it
    """

    directions: np.ndarray
    weights: np.ndarray
    pixel_area: u.Quantity | None = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        if len(self.directions) != len(self.weights):
            raise ValueError(
                f"{len(self.directions)} directions against "
                f"{len(self.weights)} weights"
            )

    def __len__(self):
        return len(self.directions)

    def __repr__(self):
        area = ("" if self.pixel_area is None
                else f", {self.area.to_value(u.deg**2):.0f} deg2")
        label = self.meta.get("label", "region")
        return f"SkyRegion({label}, {len(self)} directions{area})"

    # -- construction ------------------------------------------------------

    @classmethod
    def from_lonlat(cls, lon, lat, weights, pixel_area=None, meta=None):
        """
        A region from longitudes and latitudes on a sphere.

        Parameters
        ----------
        lon, lat: `astropy.Quantity`
            angles, one per direction
        weights: array-like
        pixel_area: `astropy.Quantity`, optional
        meta: dict, optional

        Returns
        -------
        `SkyRegion`
        """
        lon = u.Quantity(lon).to_value(u.rad)
        lat = u.Quantity(lat).to_value(u.rad)
        directions = np.column_stack([np.cos(lat) * np.cos(lon),
                                      np.cos(lat) * np.sin(lon),
                                      np.sin(lat)])
        return cls(directions, np.asarray(weights, dtype=float),
                   pixel_area, dict(meta or {}))

    @classmethod
    def from_coords(cls, coords, weights, observation=None, pixel_area=None,
                    meta=None):
        """
        A region from sky coordinates, optionally put over the array.

        Parameters
        ----------
        coords: `astropy.coordinates.SkyCoord`
        weights: array-like
        observation: `divtel.observation.Observation`, optional
            place and time. Given one, the directions come out in the array's
            horizontal frame, which is where the region really was that night.
            Without one they stay in the coordinates' own frame, ready for
            `place` to set them down wherever the question wants them.
        pixel_area: `astropy.Quantity`, optional
        meta: dict, optional

        Returns
        -------
        `SkyRegion`
        """
        if observation is not None:
            local = coords.transform_to(observation.altaz)
            directions = np.column_stack(
                alt_az_to_vector(local.alt, local.az))
            return cls(directions, np.asarray(weights, dtype=float),
                       pixel_area, dict(meta or {}))

        cartesian = coords.cartesian
        directions = np.column_stack([cartesian.x.value, cartesian.y.value,
                                      cartesian.z.value])
        return cls(directions, np.asarray(weights, dtype=float),
                   pixel_area, dict(meta or {}))

    @classmethod
    def from_table(cls, path):
        """
        A region read back from an ECSV table of ``ra``, ``dec`` and ``prob``.

        The format `divtel.data.gw170817` ships its precomputed credible regions
        in, so a study -- or a notebook running in a browser, where no HEALPix
        reader is available -- can start from a real localization without the
        sky map it was cut from.

        Parameters
        ----------
        path: str or `pathlib.Path`

        Returns
        -------
        `SkyRegion`
            in the table's own celestial frame; call `place` to set it down
        """
        from astropy.table import Table

        path = Path(path)
        if path.suffix == ".gz":
            # astropy identifies a format from the extension, and ".ecsv.gz" is
            # not one it knows; say so, and hand it a decompressed stream.
            with gzip.open(path, "rt") as handle:
                # The ascii reader wants lines, not a file object.
                table = Table.read(handle.read().splitlines(),
                                   format="ascii.ecsv")
        else:
            table = Table.read(path)
        meta = dict(table.meta)
        pixel_area = meta.pop("pixel_area_deg2", None)
        return cls.from_lonlat(
            u.Quantity(table["ra"]), u.Quantity(table["dec"]), table["prob"],
            pixel_area * u.deg**2 if pixel_area is not None else None, meta)

    def to_table(self, path=None):
        """
        The region as an ECSV table, the format `from_table` reads.

        Parameters
        ----------
        path: str or `pathlib.Path`, optional
            write it here as well as returning it

        Returns
        -------
        `astropy.table.Table`
        """
        from astropy.table import Table

        lat = np.arcsin(np.clip(self.directions[:, 2], -1.0, 1.0))
        lon = np.arctan2(self.directions[:, 1], self.directions[:, 0]) % (2 * np.pi)

        table = Table({"ra": (lon * u.rad).to(u.deg),
                       "dec": (lat * u.rad).to(u.deg),
                       "prob": self.weights})
        # Written to the precision the numbers are worth rather than to the
        # precision a float64 prints at: five decimals of a degree is a tenth of
        # an arcsecond, against pixels a fifth of a degree across. It halves the
        # file, and these files are carried into a browser.
        table["ra"].format = table["dec"].format = "%.5f"
        table["prob"].format = "%.6e"
        table.meta.update(self.meta)
        if self.pixel_area is not None:
            table.meta["pixel_area_deg2"] = float(
                self.pixel_area.to_value(u.deg**2))

        if path is not None:
            path = Path(path)
            if path.suffix == ".gz":
                with gzip.open(path, "wt") as handle:
                    table.write(handle, format="ascii.ecsv")
            else:
                table.write(path, format="ascii.ecsv", overwrite=True)
        return table

    # -- placement ---------------------------------------------------------

    def placement_rotation(self, alt, az):
        """
        The rotation that would carry this region's centroid to a pointing.

        Separate from `place` so that several regions cut from the same map can
        be moved together. Drawing a picture needs exactly that: the 50, 90 and
        99 per cent bands have different centroids, and rotating each by its own
        would slide them apart instead of nesting them.

        Parameters
        ----------
        alt, az: `astropy.Quantity`

        Returns
        -------
        `numpy.ndarray`
            a 3 by 3 rotation matrix
        """
        return rotation_between(self.centroid,
                                alt_az_to_vector(alt, az).astype(float))

    def rotate(self, rotation):
        """
        The same region turned by a rotation matrix.

        Parameters
        ----------
        rotation: `numpy.ndarray`
            3 by 3, as `placement_rotation` returns

        Returns
        -------
        `SkyRegion`
        """
        return replace(self, directions=self.directions @ rotation.T)

    def place(self, alt, az):
        """
        The region set down at a chosen place over the array, shape intact.

        For a strategy question a localization supplies a *shape*: how elongated
        it is, how its probability is distributed along the arc. Where it
        happened to be on some night is an accident nobody is planning around.
        This rotates the region rigidly so its centroid sits at the given
        altitude and azimuth, leaving every internal angle unchanged.

        Controlling the elevation is not cosmetic. An array is foreshortened
        when it looks near the horizon, so the same divergence produces a
        different spread and a different multiplicity at 20 degrees than at 70.
        Comparing regions evaluated at whatever elevation each happened to fall
        at would mix that in with the effect being measured.

        Parameters
        ----------
        alt, az: `astropy.Quantity`
            where to put the centroid

        Returns
        -------
        `SkyRegion`
            in the array's ground frame

        Notes
        -----
        The rotation is the shortest one carrying the centroid to the target, so
        the region is not spun about its own axis beyond what that implies. For
        an elongated region the position angle it lands at is therefore
        arbitrary, which is honest, since the position angle of a future alert
        is too.
        """
        return self.rotate(self.placement_rotation(alt, az))

    def visible_part(self, alt_min=0 * u.deg):
        """
        The part of a placed region a telescope could actually point at.

        `place` rotates a localization so its *centroid* sits at a chosen
        altitude, which says nothing about its ends. GW170817's first alert map
        is an arc 124 degrees long; put its centroid at 60 degrees elevation and
        a sixth of its probability is below the horizon. An optimiser handed
        that region will cheerfully aim telescopes underground and report the
        coverage, and every strategy compared against it inherits the fiction.

        Cutting the region down first is the fix, and it has to happen before
        the strategies rather than after: what is left is a different shape, not
        just a smaller one, and the right pointing for it is different too.

        Parameters
        ----------
        alt_min: `astropy.Quantity`
            lowest altitude a telescope may point at. Zero is the horizon and
            pure geometry; a real observation stops well above it, both because
            the drives do and because a shower seen through fifteen atmospheres
            is not worth pointing at.

        Returns
        -------
        region: `SkyRegion`
            the directions above `alt_min`, with their weights unchanged -- so
            they still sum to the whole region's weight, and coverage measured
            against them stays a fraction of the whole
        dropped: float
            the weight cut away, as a fraction of the region
        """
        keep = self.directions[:, 2] >= np.sin(alt_min.to_value(u.rad))
        dropped = float(self.weights[~keep].sum() / self.weights.sum())
        return self[keep], dropped

    def __getitem__(self, selection):
        """The region restricted to a subset of its directions."""
        return replace(self, directions=self.directions[selection],
                       weights=self.weights[selection])

    # -- size --------------------------------------------------------------

    @property
    def area(self):
        """
        Solid angle of the region.

        Returns
        -------
        `astropy.Quantity`
            in deg**2

        Raises
        ------
        ValueError
            if the region carries no pixel area to sum
        """
        if self.pixel_area is None:
            raise ValueError(
                "this region has no pixel area, so its solid angle is unknown; "
                "build it from a pixelized map, or pass pixel_area"
            )
        return len(self) * self.pixel_area.to(u.deg**2)

    @property
    def centroid(self):
        """
        The weight-weighted centre of the region, as a unit vector.

        The mean of the direction vectors, renormalized -- the direction the
        region sits in, treated as a whole. For an arc this lands off the arc
        itself, inside the sphere's chord, which is correct: it is a centre of
        direction, not a point of the region.
        """
        mean = np.average(self.directions, axis=0, weights=self.weights)
        return mean / np.linalg.norm(mean)

    @property
    def equivalent_radius(self):
        """
        The radius the region would have if it were a disc.

        The naive size, quoted for contrast with `containment_radius`.

        Returns
        -------
        `astropy.Quantity`
            in degrees
        """
        return np.sqrt(self.area / np.pi).to(
            u.deg, equivalencies=u.dimensionless_angles())

    def containment_radius(self, iterations=2000, axis=False):
        """
        The smallest cone that contains the whole region.

        This is the number an array has to match: point the array at the cone's
        axis, and every telescope has to reach `containment_radius` off-axis for
        the array to see the whole region at once. Twice it is the angular width
        a divergent configuration must span.

        Found by moving a trial axis repeatedly a shrinking step towards
        whichever direction is currently farthest from it -- the standard
        incremental construction for a minimum enclosing ball, run on the unit
        vectors. It converges from above onto the true minimum, and it needs no
        derivatives and no scipy. The weighted centroid is a good enough
        starting point that a couple of thousand steps settle it to well under a
        hundredth of a degree.

        Parameters
        ----------
        iterations: int
            steps of the search; more is tighter and slower
        axis: bool
            also return the cone's axis, as a unit vector

        Returns
        -------
        `astropy.Quantity`, or (`astropy.Quantity`, `numpy.ndarray`)
            the cone's angular radius, in degrees

        Notes
        -----
        A region can be more than a hemisphere across, and then no cone smaller
        than the whole sky contains it. The radius simply comes out above 90
        degrees, and the caller should read that as "no single pointing reaches
        this".
        """
        centre = np.average(self.directions, axis=0)
        centre /= np.linalg.norm(centre)

        for step in range(1, iterations + 1):
            farthest = self.directions[np.argmin(self.directions @ centre)]
            centre = centre + (farthest - centre) / (step + 1)
            centre /= np.linalg.norm(centre)

        radius = np.arccos(
            np.clip(np.min(self.directions @ centre), -1.0, 1.0)) * u.rad

        return (radius.to(u.deg), centre) if axis else radius.to(u.deg)

    @property
    def max_separation(self):
        """
        The greatest angular distance between two directions of the region.

        The region's length, as opposed to its width. For the arc-shaped maps an
        early two-detector alert produces this runs several times the equivalent
        radius, and that gap is the whole reason non-parallel pointing is worth
        considering.

        Returns
        -------
        `astropy.Quantity`
            in degrees

        Notes
        -----
        Exact up to twenty thousand directions; above that an evenly-spaced
        subsample stands in, which can only underestimate, and by less than the
        spacing between directions.
        """
        vectors = self.directions
        if len(vectors) > _MAX_EXACT_PIXELS:
            stride = int(np.ceil(len(vectors) / _MAX_EXACT_PIXELS))
            vectors = vectors[::stride]

        # Chunked so the pairwise dot products never need a matrix of the whole
        # region against itself, which for a large region would not fit.
        smallest = 1.0
        for start in range(0, len(vectors), 1024):
            smallest = min(smallest,
                           float((vectors[start:start + 1024] @ vectors.T).min()))

        return (np.arccos(np.clip(smallest, -1.0, 1.0)) * u.rad).to(u.deg)

    def describe_extent(self):
        """
        Every size measure at once.

        Returns
        -------
        dict
            ``area``, ``equivalent_radius``, ``containment_radius``,
            ``max_separation`` and ``elongation`` -- the last being the
            containment radius in units of the equivalent radius, so 1 is a disc
            and anything much above it is a stretched region that costs more to
            cover than its area suggests.
        """
        equivalent = self.equivalent_radius
        containment, axis = self.containment_radius(axis=True)
        return {
            "area": self.area,
            "weight": float(self.weights.sum()),
            "equivalent_radius": equivalent,
            "containment_radius": containment,
            "containment_axis": axis,
            "max_separation": self.max_separation,
            "elongation": float(containment / equivalent),
        }

    # -- what an array sees ------------------------------------------------

    def multiplicity(self, array):
        """
        How many of the array's cameras cover each direction.

        A direction is seen by a telescope when it lies within that telescope's
        field-of-view radius of where the telescope points. This is the quantity
        every strategy in `divtel.strategy` is spending: a shower seen once
        cannot be reconstructed stereoscopically, so multiplicity and not
        coverage is what an array actually delivers at a direction.

        Parameters
        ----------
        array: `divtel.telescope.Array`
            already pointed

        Returns
        -------
        `numpy.ndarray` of int
            one entry per direction
        """
        cos_radii = np.cos([t.fov_radius.to_value(u.rad)
                            for t in array.telescopes])
        # One matrix product rather than a pass per telescope: the sweeps in
        # `divtel.strategy` evaluate this thousands of times.
        return ((self.directions @ array.pointing_vectors.T)
                >= cos_radii).sum(axis=1)

    def covered(self, array, m_cut=2):
        """
        The weight of the region seen by at least `m_cut` telescopes.

        The `m_cut` is the whole reason spreading an array is a trade rather
        than a free win. A shower seen by one telescope cannot be reconstructed
        stereoscopically, so ``m_cut=2`` is what a real analysis needs, and it
        is exactly what spreading eats into.

        Parameters
        ----------
        array: `divtel.telescope.Array`
            already pointed
        m_cut: int

        Returns
        -------
        float
            as a fraction of the region's total weight
        """
        weights = self.weights
        return float(weights[self.multiplicity(array) >= m_cut].sum()
                     / weights.sum())

    def visible_fraction(self, observation, min_altitude=0 * u.deg):
        """
        How much of the region's weight is above the horizon at a given time.

        For a region in celestial coordinates, not a placed one. A localization
        tens of degrees across sets over an hour or more, so for part of the
        night the array can only reach part of it however it points. This is the
        ceiling on what any pointing can cover.

        Parameters
        ----------
        observation: `divtel.observation.Observation`
        min_altitude: `astropy.Quantity`

        Returns
        -------
        float
            between 0 and 1
        """
        lat = np.arcsin(np.clip(self.directions[:, 2], -1.0, 1.0)) * u.rad
        lon = np.arctan2(self.directions[:, 1], self.directions[:, 0]) * u.rad
        coords = SkyCoord(lon, lat, frame="icrs")

        altitude = coords.transform_to(observation.altaz).alt
        up = altitude >= min_altitude
        return float(self.weights[up].sum() / self.weights.sum())

    @property
    def coord(self):
        """
        The region's directions as an `astropy.coordinates.SkyCoord`, in ICRS.

        Only meaningful for a region still in celestial coordinates; a placed
        region's vectors are directions over the array, and reading them as
        right ascension and declination is nonsense.
        """
        lat = np.arcsin(np.clip(self.directions[:, 2], -1.0, 1.0)) * u.rad
        lon = np.arctan2(self.directions[:, 1], self.directions[:, 0]) * u.rad
        return SkyCoord(lon, lat, frame="icrs")

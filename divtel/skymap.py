"""
Read a gravitational-wave localization and cut credible regions out of it.

This module is the optional half of divtel. It needs a HEALPix reader, which
nothing else here does, so it lives behind an extra::

    pip install divtel[skymap]

What it produces is a `divtel.region.SkyRegion` -- unit vectors and weights --
and from there `divtel.strategy` and everything else work with no knowledge of
sky maps at all. Precomputed regions for GW170817 ship in ``divtel/data``, so a
study that only wants those needs neither this module nor the extra.

A GW sky map is a probability distribution over the whole sky, stored as a
HEALPix array: the sky is split into equal-area pixels and each carries the
probability that the source is in it. Nothing in it is a position -- the source
is somewhere, and the map says where it is likely to be. What an observer needs
instead is a *region*: the smallest patch of sky holding, say, 90% of the
probability, so that pointing at it has a 90% chance of pointing at the source.

That patch comes from the greedy credible-level construction in
`credible_levels`, and everything downstream -- the angular size in
`divtel.region.SkyRegion`, the coverage in `divtel.strategy` -- works from it.

The maps this reads are the ones LIGO/Virgo actually circulated in 2017, listed
in `SKYMAPS`. They are flat HEALPix arrays in the NESTED ordering; the
multi-order (UNIQ) format that later alerts use came afterwards, and `load`
rejects it with a message rather than reading it wrongly.
"""

from __future__ import annotations

import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

try:
    from astropy_healpix import HEALPix, nside_to_pixel_area
except ImportError as error:  # pragma: no cover - exercised by the extra
    raise ImportError(
        "divtel.skymap reads HEALPix sky maps and needs astropy-healpix, which "
        "the rest of divtel does not; install it with `pip install divtel[skymap]`"
    ) from error

from .region import SkyRegion

__all__ = [
    "SKYMAPS",
    "SkyMapInfo",
    "SkyMap",
    "fetch",
    "load",
    "credible_levels",
    "credible_mask",
    "credible_area",
    "split_regions",
    "pixel_coords",
    "region",
]

# The base of the LIGO document holding the GW170817 localizations, LIGO-G1701985
# ("GW170817 sky localization"). The 2017 alerts pointed at GraceDB, which now
# refuses anonymous access; this document is the public copy of the same files.
_DCC = "https://dcc.ligo.org/public/0146/G1701985/001/"


@dataclass(frozen=True)
class SkyMapInfo:
    """
    One published localization of an event, and when observers first had it.

    The time matters as much as the map. A follow-up is not planned against the
    best localization ever produced but against the best one available when the
    telescope could move, and for GW170817 those are different maps.

    Attributes
    ----------
    url: str
        where to download it
    filename: str
        what to cache it as
    label: str
        short human-readable name, used in tables and plot legends
    gcn: str
        the GCN circular that carried it
    t_available: str
        UTC time that circular went out, ISO format
    detectors: str
        the detector network behind it
    description: str
    """

    url: str
    filename: str
    label: str
    gcn: str
    t_available: str
    detectors: str
    description: str


# The 2017 alert sequence, in the order observers received it, plus the final
# catalogue map for the follow-on study. Circular times are from the GCN archive
# at https://gcn.gsfc.nasa.gov/other/170817A.gcn3.
SKYMAPS = {
    "bayestar_hl": SkyMapInfo(
        url=_DCC + "bayestar_no_virgo.fits.gz",
        filename="bayestar_no_virgo.fits.gz",
        label="BAYESTAR H+L",
        gcn="GCN 21509",
        t_available="2017-08-17T13:21:00",
        detectors="H1,L1",
        description="First rapid localization, before Virgo data were folded in.",
    ),
    "bayestar_hlv": SkyMapInfo(
        url=_DCC + "bayestar.fits.gz",
        filename="bayestar.fits.gz",
        label="BAYESTAR H+L+V",
        gcn="GCN 21513",
        t_available="2017-08-17T17:54:00",
        detectors="H1,L1,V1",
        description="Rapid localization with all three detectors; the map that "
                    "sent observers to NGC 4993.",
    ),
    "lalinference_prelim": SkyMapInfo(
        url=_DCC + "preliminary-LALInference.fits.gz",
        filename="preliminary-LALInference.fits.gz",
        label="LALInference (preliminary)",
        gcn="GCN 21527",
        t_available="2017-08-18T05:00:00",
        detectors="H1,L1,V1",
        description="First full parameter-estimation localization.",
    ),
    "lalinference_v2": SkyMapInfo(
        url=_DCC + "LALInference_v2.fits.gz",
        filename="LALInference_v2.fits.gz",
        label="LALInference v2",
        gcn="GCN 21983",
        t_available="2017-08-30T00:00:00",
        detectors="H1,L1,V1",
        description="Updated parameter estimation, two weeks after the merger.",
    ),
    "gwtc1": SkyMapInfo(
        url="https://dcc.ligo.org/public/0157/P1800381/007/GW170817_skymap.fits.gz",
        filename="GW170817_skymap.fits.gz",
        label="GWTC-1 final",
        gcn="GWTC-1",
        t_available="2018-11-30T00:00:00",
        detectors="H1,L1,V1",
        description="Final catalogue localization. Not part of the alert "
                    "sequence -- kept here for the follow-on study.",
    ),
}


@dataclass
class SkyMap:
    """
    A HEALPix probability sky map.

    Attributes
    ----------
    prob: `numpy.ndarray`
        probability per pixel, summing to one
    nside: int
        HEALPix resolution
    order: str
        ``"nested"`` or ``"ring"``
    info: `SkyMapInfo` or None
        which published map this is, when it came from `SKYMAPS`
    meta: dict
        the FITS header, for provenance
    """

    prob: np.ndarray
    nside: int
    order: str = "nested"
    info: SkyMapInfo | None = None
    meta: dict = field(default_factory=dict)

    @property
    def label(self):
        """Short name for tables and legends."""
        return self.info.label if self.info else "sky map"

    @property
    def healpix(self):
        """The `astropy_healpix.HEALPix` describing this pixelization."""
        return HEALPix(nside=self.nside, order=self.order, frame="icrs")

    @property
    def pixel_area(self):
        """
        Solid angle of one pixel

        Returns
        -------
        `astropy.Quantity`
            in deg**2
        """
        return nside_to_pixel_area(self.nside).to(u.deg**2)

    def degrade(self, nside):
        """
        The same map at a coarser resolution.

        A nside-2048 map is fifty million pixels, and the geometry downstream
        touches every one of them in the credible region. Coarsening first keeps
        that affordable. In the NESTED ordering the four children of a pixel are
        adjacent in the array, so degrading is a reshape and a sum -- and since
        probability is summed rather than averaged, the total is preserved
        exactly and credible areas move by less than one coarse pixel.

        Parameters
        ----------
        nside: int
            target resolution, a power of two no larger than the current one

        Returns
        -------
        `SkyMap`

        Raises
        ------
        ValueError
            if the target is finer than the map, or not a power of two
        """
        if nside == self.nside:
            return self
        if nside > self.nside:
            raise ValueError(
                f"cannot degrade an nside-{self.nside} map to nside {nside}; "
                "degrading only makes a map coarser"
            )
        if nside < 1 or nside & (nside - 1):
            raise ValueError(f"nside must be a power of two, got {nside}")
        if self.order != "nested":
            raise ValueError(
                "degrading needs the NESTED ordering, where a pixel's children "
                f"are adjacent in the array; this map is {self.order.upper()}"
            )

        factor = (self.nside // nside) ** 2
        prob = self.prob.reshape(-1, factor).sum(axis=1)
        return SkyMap(prob=prob, nside=nside, order=self.order,
                      info=self.info, meta=self.meta)


def fetch(name, cache_dir, download=True):
    """
    The file for a published sky map, downloading it once.

    Parameters
    ----------
    name: str
        a key of `SKYMAPS`
    cache_dir: str or `pathlib.Path`
        directory to keep downloads in; created if missing
    download: bool
        fetch a missing file. With False a missing file is an error, which is
        what you want when re-running offline and a silent download would hang.

    Returns
    -------
    `pathlib.Path`

    Raises
    ------
    KeyError
        if the name is not one of the published maps
    FileNotFoundError
        if the file is missing and downloading is off
    """
    if name not in SKYMAPS:
        raise KeyError(
            f"unknown sky map {name!r}; known maps are {', '.join(SKYMAPS)}"
        )

    info = SKYMAPS[name]
    cache_dir = Path(cache_dir)
    path = cache_dir / info.filename

    if path.exists():
        return path

    if not download:
        raise FileNotFoundError(
            f"{path} is not cached and downloading is off; drop --no-download "
            f"to fetch it from {info.url}"
        )

    cache_dir.mkdir(parents=True, exist_ok=True)
    # Download beside the target and rename, so an interrupted download cannot
    # leave a truncated file that later runs would happily read as a sky map.
    partial = path.with_suffix(path.suffix + ".part")
    urllib.request.urlretrieve(info.url, partial)
    partial.rename(path)

    return path


def load(path, info=None, max_nside=None):
    """
    Read a flat HEALPix sky map from a FITS file.

    Parameters
    ----------
    path: str or `pathlib.Path`
    info: `SkyMapInfo`, optional
        which published map this is
    max_nside: int, optional
        coarsen the map to at most this resolution, via `SkyMap.degrade`

    Returns
    -------
    `SkyMap`

    Raises
    ------
    ValueError
        if the file is a multi-order (UNIQ) map, or has no probability column
    """
    table = Table.read(path)

    if "UNIQ" in table.colnames:
        raise ValueError(
            f"{path} is a multi-order (UNIQ) sky map, which this reader does not "
            "handle; the 2017 GW170817 alerts are all flat maps"
        )

    if "PROB" not in table.colnames:
        raise ValueError(
            f"{path} has no PROB column; its columns are "
            f"{', '.join(table.colnames)}"
        )

    prob = np.asarray(table["PROB"], dtype=float)
    order = str(table.meta.get("ORDERING", "NESTED")).strip().lower()
    if order not in ("nested", "ring"):
        raise ValueError(f"{path} has an unrecognised ORDERING {order!r}")

    nside = int(table.meta.get("NSIDE", np.sqrt(len(prob) / 12)))

    skymap = SkyMap(prob=prob, nside=nside, order=order, info=info,
                    meta=dict(table.meta))

    if max_nside is not None and nside > max_nside:
        skymap = skymap.degrade(max_nside)

    return skymap


def credible_levels(prob):
    """
    The credible level at which each pixel joins the region.

    Build the region greedily: take the most probable pixel, then the next, and
    so on. A pixel's credible level is the total probability accumulated by the
    time it is taken, so the 90% region is every pixel with a level below 0.9.
    This is the standard construction, and it gives the *smallest* region
    holding a given probability -- any other region of the same area holds less.

    Parameters
    ----------
    prob: `numpy.ndarray`
        probability per pixel

    Returns
    -------
    `numpy.ndarray`
        credible level per pixel, in the same order, between 0 and 1
    """
    order = np.argsort(prob)[::-1]
    levels = np.empty(len(prob))
    levels[order] = np.cumsum(prob[order])
    return levels


def credible_mask(skymap, level=0.9):
    """
    Which pixels are in the smallest region holding this much probability.

    Parameters
    ----------
    skymap: `SkyMap`
    level: float
        probability to contain, e.g. 0.9 for the 90% credible region

    Returns
    -------
    `numpy.ndarray` of bool
        one entry per pixel
    """
    if not 0 < level <= 1:
        raise ValueError(f"level must be in (0, 1], got {level}")
    return credible_levels(skymap.prob) <= level


def credible_area(skymap, level=0.9):
    """
    Area of the credible region.

    Parameters
    ----------
    skymap: `SkyMap`
    level: float

    Returns
    -------
    `astropy.Quantity`
        in deg**2
    """
    return credible_mask(skymap, level).sum() * skymap.pixel_area


def split_regions(skymap, level=0.9):
    """
    A credible region cut into its separate pieces.

    A credible region need not be one patch of sky. With two detectors the
    timing leaves an ambiguity that often puts a second lobe somewhere else
    entirely, and GW170817's first map is a case in point: a lobe holding 61% of
    the probability, and another holding 29% eighty degrees away.

    That matters for planning, because no array covers both at once and treating
    them as one region asks it to. An observer picks a lobe. This returns them
    separately so a study can do the same.

    Pieces are found by walking the region through HEALPix pixel adjacency, so
    two pixels are in the same piece exactly when the region connects them.

    Parameters
    ----------
    skymap: `SkyMap`
    level: float
        credible level to cut at

    Returns
    -------
    list of `numpy.ndarray` of bool
        one pixel mask per piece, largest probability first. Their union is
        `credible_mask(skymap, level)`.

    Examples
    --------
    >>> pieces = split_regions(skymap, 0.9)
    >>> [round(float(skymap.prob[m].sum()), 3) for m in pieces]
    [0.614, 0.286]
    """
    mask = credible_mask(skymap, level)
    healpix = skymap.healpix
    inside = np.flatnonzero(mask)
    membership = {int(pixel): index for index, pixel in enumerate(inside)}

    # Union-find over the region's pixels, joining each to any neighbour that is
    # also in the region. One pass, and no recursion to overflow on a region
    # that wraps most of the way round the sky.
    parent = np.arange(len(inside))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    neighbours = healpix.neighbours(inside)
    for position, pixel_neighbours in enumerate(neighbours.T):
        for neighbour in pixel_neighbours:
            other = membership.get(int(neighbour))
            if other is None:
                continue
            a, b = find(position), find(other)
            if a != b:
                parent[a] = b

    roots = np.array([find(i) for i in range(len(inside))])

    pieces = []
    for root in np.unique(roots):
        piece = np.zeros(len(skymap.prob), dtype=bool)
        piece[inside[roots == root]] = True
        pieces.append(piece)

    pieces.sort(key=lambda piece: -skymap.prob[piece].sum())
    return pieces


def pixel_coords(skymap, mask=None):
    """
    Where the pixels are on the sky.

    Parameters
    ----------
    skymap: `SkyMap`
    mask: `numpy.ndarray` of bool, optional
        take only these pixels; all of them by default

    Returns
    -------
    `astropy.coordinates.SkyCoord`
        in ICRS, one entry per selected pixel
    """
    index = np.arange(len(skymap.prob)) if mask is None else np.flatnonzero(mask)
    lon, lat = skymap.healpix.healpix_to_lonlat(index)
    return SkyCoord(lon, lat, frame="icrs")


def region(skymap, level=0.9, mask=None):
    """
    A credible region of a sky map, as a `divtel.region.SkyRegion`.

    The bridge out of this module. Everything downstream -- the sizes in
    `divtel.region.SkyRegion`, the strategies in `divtel.strategy` -- works on
    unit vectors and weights, and past this call nothing knows the region came
    from a HEALPix map.

    The vectors come out in the map's own celestial frame, so the region is
    still where the sky put it. Use `divtel.region.SkyRegion.place` to set it
    down over an array at a chosen elevation, or pass an
    `divtel.observation.Observation` to
    `divtel.region.SkyRegion.from_coords` instead if the question is about a
    particular night.

    Parameters
    ----------
    skymap: `SkyMap`
    level: float
        credible level to cut at, e.g. 0.9
    mask: `numpy.ndarray` of bool, optional
        take these pixels instead of the credible region, so one piece of a
        region `split_regions` has separated can be carried on alone

    Returns
    -------
    `divtel.region.SkyRegion`
    """
    if mask is None:
        mask = credible_mask(skymap, level)

    info = skymap.info
    meta = {"level": float(level)}
    if info is not None:
        meta.update({"label": info.label, "gcn": info.gcn,
                     "t_available": info.t_available,
                     "detectors": info.detectors})

    return SkyRegion.from_coords(pixel_coords(skymap, mask), skymap.prob[mask],
                                 pixel_area=skymap.pixel_area, meta=meta)

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from . import pointing


class Telescope:
    """
    Describe a generic telescope

    Parameters
    ----------
    x: `astropy.quantity`
        x ground position
    y: `astropy.quantity`
        y ground position
    z: `astropy.quantity`
        z ground position
    focal: `astropy.quantity`
    camera_radius: `astropy.quantity`
    tel_id: int, optional
        telescope identifier. CTAO numbers telescopes by type -- LSTs 1 to 4,
        MSTs 5 to 14, telescopes added to try out a configuration from 15 on --
        so the id carries meaning and is worth setting explicitly; `Array.group_by`
        selects on it. Left out, ids come from a counter shared by every
        `Telescope` ever built, which makes them unique but not meaningful.
    """

    _id = 0

    def __init__(self, x, y, z, focal, camera_radius, tel_id=None):

        self.x = x.to(u.m)
        self.y = y.to(u.m)
        self.z = z.to(u.m)
        self.focal = focal.to(u.m)
        self.camera_radius = camera_radius.to(u.m)
        self.alt = u.Quantity(0, u.rad)
        self.az = u.Quantity(0, u.rad)
        if tel_id is None:
            Telescope._id += 1
            tel_id = Telescope._id
        self.id = tel_id

    def point_to_altaz(self, alt, az):
        """
        Point the telescope in a given direction

        Parameters
        ----------
        alt: `astropy.Quantity`
            altitude, from the ground towards z
        az: `astropy.Quantity`
            azimuth, clock-wise from x towards y
        """
        self.alt = alt.to(u.rad)
        self.az = az.to(u.rad)

    @property
    def zenith(self):
        """
        Zenith angle of the current pointing

        Returns
        -------
        `astropy.Quantity`
            angle from the vertical, in radians; the complement of `alt`
        """
        return np.pi/2.*u.rad - self.alt

    @property
    def fov(self):
        """
        Area of the field of view in rad**2
        """
        # camera_radius / focal is dimensionless, so it has to be stripped of
        # its unit before being declared an angle -- astropy will not convert
        # dimensionless to rad**2 on its own.
        half_angle = (self.camera_radius / self.focal).to_value(u.dimensionless_unscaled)
        return u.Quantity(np.pi * half_angle**2, u.rad**2)

    @property
    def fov_radius(self):
        """
        Angular radius of the field of view

        Returns
        -------
        `astropy.Quantity`
            half-angle subtended by the camera, in radians
        """
        return u.Quantity(np.arctan2(self.camera_radius, self.focal), u.rad)

    @property
    def position(self):
        """
        Ground position of the telescope

        Returns
        -------
        `astropy.Quantity`
            [x, y, z] in metres
        """
        return np.array([self.x.to_value(u.m),
                         self.y.to_value(u.m),
                         self.z.to_value(u.m)]) * u.m

    def point_to_object(self, target):
        """
        Point to object.

        Parameters
        ----------
        target: `astropy.Quantity` or numpy.array([x, y, z])
            position of the object; a plain array is read as metres
        """
        target = u.Quantity(target, u.m)
        GT = np.sqrt(((self.position - target) ** 2).sum())
        alt_tel = np.arcsin((target[2] - self.z) / GT)
        # Az runs clock-wise from X towards Y (see divtel.pointing), so the
        # Y offset enters negated: az = arctan2(-dy, dx).
        az_tel = np.arctan2(self.y - target[1], target[0] - self.x)
        self.point_to_altaz(alt_tel.to(u.rad), az_tel.to(u.rad))

    @property
    def pointing_vector(self):
        """
        Current pointing direction as a unit vector

        Returns
        -------
        `numpy.array`
            [x, y, z], unit length, in the ground frame
        """
        # return pointing.alt_az_to_vector(self.alt, self.az)
        return np.array([np.cos(self.alt.to(u.rad))*np.cos(self.az.to(u.rad)),
                         -np.cos(self.alt.to(u.rad))*np.sin(self.az.to(u.rad)),
                         np.sin(self.alt.to(u.rad))])


class Array:
    """
    Describe a telescope array

    Parameters
    ----------
    telescope_list: [Telescope]
        list of telescopes forming the array
    """

    def __init__(self, telescope_list):
        self.telescopes = telescope_list
        # The divergence and the direction diverged from are not recoverable
        # from the telescopes afterwards -- G lies behind the array, so the
        # mean of the telescope pointings is not the requested mean pointing.
        # `export_cfg` needs both, so record them when they are handed to us.
        self._div = 0
        self._alt_mean = None
        self._az_mean = None

    @property
    def positions_array(self):
        """
        All telescopes positions

        Returns
        -------
        positions_array: `astropy.Quantity`
            2D array of [x, y, z] in metres
        """
        return u.Quantity([tel.position for tel in self.telescopes])

    @property
    def pointing_vectors(self):
        """
        all telescopes pointing vectors as an array

        Returns
        -------
        pointing_vectors: `numpy.array`
            2D numpy array
        """
        return np.array([tel.pointing_vector for tel in self.telescopes])

    @property
    def barycenter(self):
        """
        Barycenter of the array

        Returns
        -------
        `astropy.Quantity`
            [x, y, z] in metres
        """
        return self.positions_array.mean(axis=0)

    @property
    def pointing_altaz(self):
        """
        All telescopes pointing directions (alt, az) as an array

        Returns
        -------
        pointing_directions: `astropy.Quantity`
            2D array: [[alt1,az1], [alt2, az2], ...]
        """
        return np.array([np.array([tel.alt.to_value(u.rad), tel.az.to_value(u.rad)]) for tel in self.telescopes])*u.rad

    @property
    def div(self):
        """
        Divergence parameter last passed to `divergent_pointing`

        Returns
        -------
        float
        """
        return self._div

    @property
    def mean_pointing(self):
        """
        Direction the array was last asked to diverge from

        Returns
        -------
        (alt, az): tuple of `astropy.Quantity`, in degrees

        Raises
        ------
        ValueError
            if the array has never been pointed
        """
        if self._alt_mean is None:
            raise ValueError("the array has no pointing yet; "
                             "call divergent_pointing first")
        return self._alt_mean, self._az_mean

    def group_by(self, groups):
        """
        Split the array into sub-arrays.

        A real array is not one instrument but several sharing a site: on the
        CTA La Palma layout four LSTs sit inside fifteen MSTs, with different
        cameras, different fields of view, and their own barycenters. Grouping
        lets each be looked at on its own.

        Parameters
        ----------
        groups: dict or str
            either a mapping of group name to the telescope ids in it, e.g.
            ``{"LST": [1, 2, 3, 4], "MST": range(5, 20)}``, or the string
            ``"camera_radius"`` to group telescopes by the camera they carry.

        Returns
        -------
        dict of str to `Array`
            one sub-array per group, in the order the groups were given

        Raises
        ------
        ValueError
            if an id is not in the array, or a telescope is put in two groups

        Notes
        -----
        Ids are `Telescope.id`, not positions in the array. CTAO numbers
        telescopes by type -- LSTs 1 to 4, MSTs 5 to 14, telescopes added to try
        out a configuration from 15 on -- and `divtel.layout.load_array` takes
        those ids from the layout file, so they mean what the layout says.

        The sub-arrays hold the *same* `Telescope` objects as this array rather
        than copies, so pointing the whole array points every group with it.
        Each sub-array is handed the array's current pointing, so `hyper_fov`
        and `export_cfg` work on it straight away.

        Telescopes may be left out of every group; nothing requires the groups
        to cover the array.

        Examples
        --------
        >>> groups = array.group_by({"LST": range(1, 5), "MST": range(5, 20)})
        >>> groups["LST"].barycenter
        <Quantity [ 0.895 , 44.8475, 44.925 ] m>
        """
        if isinstance(groups, str):
            if groups != "camera_radius":
                raise ValueError(
                    f"group_by takes a dict of ids or 'camera_radius', got {groups!r}"
                )
            return self._group_by_camera_radius()

        by_id = {tel.id: tel for tel in self.telescopes}
        seen = {}
        grouped = {}

        for name, ids in groups.items():
            ids = list(ids)

            unknown = [tel_id for tel_id in ids if tel_id not in by_id]
            if unknown:
                raise ValueError(
                    f"group {name!r} asks for telescope id(s) "
                    f"{', '.join(map(str, unknown))}, which the array does not have; "
                    f"its ids are {', '.join(str(tel.id) for tel in self.telescopes)}"
                )

            for tel_id in ids:
                if tel_id in seen:
                    raise ValueError(
                        f"telescope {tel_id} is in both group {seen[tel_id]!r} "
                        f"and group {name!r}"
                    )
                seen[tel_id] = name

            grouped[name] = self._sub_array([by_id[tel_id] for tel_id in ids])

        return grouped

    def _group_by_camera_radius(self):
        """
        One sub-array per distinct camera radius, largest camera first.

        Telescopes of a type share a camera, so on a real layout this recovers
        the types without being told them. Groups are named after the radius,
        since that is all this knows -- pass a dict to name them yourself.
        """
        radii = sorted({tel.camera_radius for tel in self.telescopes}, reverse=True)

        return {
            f"{radius.to_value(u.m):.4g} m": self._sub_array(
                [tel for tel in self.telescopes if tel.camera_radius == radius]
            )
            for radius in radii
        }

    def _sub_array(self, telescopes):
        """An `Array` of some of these telescopes, keeping the pointing state."""
        sub = Array(telescopes)
        sub._div = self._div
        sub._alt_mean = self._alt_mean
        sub._az_mean = self._az_mean
        return sub

    def divergent_pointing(self, div, alt_mean, az_mean):
        """
        Divergent pointing given a parameter div.
        Update pointing of all telescopes of the array

        Parameters
        ----------
        div: float between 0 and 1
        alt_mean: `astropy.Quantity`
            mean alt pointing
        az_mean: `astropy.Quantity`
            mean az pointing

        Raises
        ------
        ValueError
            if div is outside [0, 1]
        """
        # Outside [0, 1] the geometry stops meaning anything: div > 1 makes
        # arcsin return nan (a RuntimeWarning at most), and div < 0 puts G in
        # front of the array rather than behind it, so every telescope ends up
        # pointing below the horizon while the caller asked for alt_mean.
        if not 0 <= div <= 1:
            raise ValueError(f"div must be between 0 and 1, got {div}")

        self._div = div
        self._alt_mean = alt_mean.to(u.deg)
        self._az_mean = az_mean.to(u.deg)

        if div == 0:
            for tel in self.telescopes:
                tel.point_to_altaz(alt_mean, az_mean)
        else:
            g_point = pointing.pointG_position(self.barycenter, div, alt_mean, az_mean)
            for tel in self.telescopes:
                alt_tel, az_tel = pointing.tel_div_pointing(tel.position, g_point)
                tel.point_to_altaz(alt_tel, az_tel)

    def hyper_fov(self, m_cut=1, rim_points=128):
        """
        Hyper field of view: the sky area covered by the array's cameras.

        Each telescope sees a disc on the sky. Pointing divergently spreads
        those discs out, trading depth of coverage for width: where several
        discs overlap a shower is seen by several telescopes and can be
        reconstructed stereoscopically, where only one covers the sky it
        cannot. The union of the discs is cut into the patches formed by their
        boundaries, and each patch is labelled with its multiplicity -- how
        many telescopes see it.

        The geometry is done on the sphere, via a Lambert azimuthal equal-area
        projection centred on the array's mean pointing. That projection
        preserves area exactly, so the patch areas are true solid angles, and
        it has no pole -- working directly in the (azimuth, altitude) plane
        instead breaks down near the zenith, where telescopes a fraction of a
        degree apart on the sky are hundreds of degrees apart in azimuth.

        Parameters
        ----------
        m_cut: int
            only count patches seen by at least this many telescopes.
            The default of 1 measures the whole covered area; 2 measures the
            part that can be reconstructed stereoscopically.
        rim_points: int
            number of points sampled around each camera's rim. The default is
            accurate to about 0.01% of a disc's area.

        Returns
        -------
        area: `astropy.Quantity`
            covered area in deg**2, counting only patches above `m_cut`
        patches: list of (`shapely.Polygon`, int)
            each patch and the number of telescopes seeing it. Coordinates are
            degrees of offset from the array's mean pointing, in the
            equal-area projection, so polygon areas are in deg**2.
        """
        from shapely.geometry import LineString, Polygon
        from shapely.ops import polygonize, unary_union

        directions = self.pointing_vectors
        radii = np.array([tel.fov_radius.to_value(u.rad) for tel in self.telescopes])

        # Centre the projection on the mean pointing. If the array points every
        # which way the mean can vanish, and any direction is as good as
        # another; the zenith is the natural choice.
        centre = directions.mean(axis=0)
        norm = np.linalg.norm(centre)
        centre = centre / norm if norm > 1e-9 else np.array([0.0, 0.0, 1.0])

        # Any two vectors completing an orthonormal frame with `centre`.
        seed = np.array([1.0, 0.0, 0.0])
        if abs(centre @ seed) > 0.9:
            seed = np.array([0.0, 1.0, 0.0])
        east = np.cross(centre, seed)
        east /= np.linalg.norm(east)
        north = np.cross(centre, east)

        def project(vectors):
            """Lambert azimuthal equal-area, in degrees from the centre."""
            along = vectors @ centre
            if np.any(along <= -1 + 1e-9):
                raise ValueError(
                    "the array points in opposite directions; no single "
                    "projection of the sky can hold it"
                )
                
            scale = np.sqrt(2.0 / (1.0 + along)) * np.degrees(1.0)
            return np.column_stack([scale * (vectors @ east),
                                    scale * (vectors @ north)])

        # Sample each camera's rim on the sphere and project it, rather than
        # drawing a circle in the projection: away from the centre the image of
        # a circle is not one, and it is the true shape we need.
        angles = np.linspace(0, 2 * np.pi, rim_points, endpoint=False)
        discs = []
        for direction, radius in zip(directions, radii):
            a = np.cross(direction, centre)
            if np.linalg.norm(a) < 1e-9:
                a = east.copy()
            a /= np.linalg.norm(a)
            b = np.cross(direction, a)

            rim = (np.cos(radius) * direction
                   + np.sin(radius) * (np.cos(angles)[:, None] * a
                                       + np.sin(angles)[:, None] * b))
            discs.append(Polygon(project(rim)))

        # Cutting the union of the *boundaries* gives the arrangement: every
        # region bounded by disc edges, each of which lies wholly inside or
        # wholly outside each disc.
        boundaries = unary_union([LineString(d.exterior.coords) for d in discs])
        patches = list(polygonize(boundaries))

        labelled = []
        for patch in patches:
            probe = patch.representative_point()
            multiplicity = sum(1 for disc in discs if disc.contains(probe))
            # A ring of discs encloses a hole, and polygonize returns that hole
            # as a patch like any other. No telescope sees it, so it is not
            # part of the field of view at all.
            if multiplicity > 0:
                labelled.append((patch, multiplicity))

        area = sum(p.area for p, m in labelled if m >= m_cut)
        return u.Quantity(area, u.deg ** 2), labelled

    def multiplicity_profile(self, patches=None):
        """
        How much sky is seen by exactly one telescope, by exactly two, and so on.

        `hyper_fov` reports how much sky the array covers; this reports how well
        it covers it. Both matter, and divergence trades one for the other.

        Parameters
        ----------
        patches: list of (`shapely.Polygon`, int), optional
            the patches from a previous `hyper_fov` call, to save computing them
            again; recomputed if not given

        Returns
        -------
        multiplicity: `numpy.ndarray` of int
            the multiplicities present, in increasing order
        area: `astropy.Quantity`
            sky seen by exactly that many telescopes, in deg**2

        Examples
        --------
        >>> array.divergent_pointing(0.05, 70 * u.deg, 180 * u.deg)
        >>> multiplicity, area = array.multiplicity_profile()
        >>> multiplicity
        array([1, 2, 3, 4, 5])
        >>> np.round(area, 1)
        <Quantity [242. , 131.3,  70.8,   8. ,   0.6] deg2>
        """
        if patches is None:
            _, patches = self.hyper_fov()

        multiplicities = sorted({multiplicity for _, multiplicity in patches})
        areas = [
            sum(patch.area for patch, m in patches if m == multiplicity)
            for multiplicity in multiplicities
        ]

        return np.array(multiplicities, dtype=int), u.Quantity(areas, u.deg ** 2)

    def multiplicity_moments(self, patches=None):
        """
        Mean and variance of the multiplicity, weighted by area.

        One number for how divergent a configuration is, and one for how evenly.
        A parallel array has every patch seen by every telescope, so the mean is
        the number of telescopes and the variance is zero; as divergence grows
        the mean falls towards one.

        Weighting is by solid angle, so a patch counts for as much sky as it
        covers -- an unweighted mean over patches would let a sliver of overlap
        count for as much as the whole field.

        Parameters
        ----------
        patches: list of (`shapely.Polygon`, int), optional
            the patches from a previous `hyper_fov` call, to save computing them
            again; recomputed if not given

        Returns
        -------
        (mean, variance): tuple of float

        Examples
        --------
        >>> array.divergent_pointing(0.05, 70 * u.deg, 180 * u.deg)
        >>> mean, variance = array.multiplicity_moments()
        >>> f"{mean:.2f} +- {np.sqrt(variance):.2f}"
        '1.66 +- 0.81'
        """
        multiplicity, area = self.multiplicity_profile(patches=patches)

        weights = area.to_value(u.deg ** 2)
        mean = np.average(multiplicity, weights=weights)
        variance = np.average((multiplicity - mean) ** 2, weights=weights)

        return float(mean), float(variance)

    def export_cfg(self, filename=None, outdir="./", tel_configs=None,
                   verbose=False):
        """
        Write the array pointing as a sim_telarray configuration file.

        sim_telarray selects a telescope by defining `TELESCOPE`, so the file
        is one `#if`/`#elif` chain: block 0 holds the array-wide defaults and
        blocks 1..N hold one telescope each.

        Angles are converted to sim_telarray's convention, which is not
        divtel's: it takes a zenith angle rather than an altitude, and its
        azimuth runs the opposite way round, so

            TELESCOPE_THETA = 90 - alt
            TELESCOPE_PHI   = (360 - az) mod 360

        Parameters
        ----------
        filename: str, optional
            defaults to 'CTA-ULTRA6-LaPalma-divX-azX-altX.cfg', built from the
            divergence and mean pointing of the array
        outdir: str or `pathlib.Path`, optional
            directory to write into; must exist
        tel_configs: [str], optional
            config file each telescope includes, in array order. Defaults
            to the La Palma layout the writer was built for: the first four
            telescopes are LSTs, the rest MSTs with NectarCam.
        verbose: bool, optional
            echo the file to stdout once written

        Returns
        -------
        `pathlib.Path`
            the file that was written

        Raises
        ------
        ValueError
            if the array has never been pointed, or `tel_configs` does not
            have one entry per telescope
        """
        from pathlib import Path

        alt_mean, az_mean = self.mean_pointing

        if tel_configs is None:
            tel_configs = ["CTA-ULTRA6-LST.cfg" if n <= 4
                           else "CTA-ULTRA6-MST-NectarCam.cfg"
                           for n in range(1, len(self.telescopes) + 1)]
        elif len(tel_configs) != len(self.telescopes):
            raise ValueError(
                f"tel_configs has {len(tel_configs)} entries but the array has "
                f"{len(self.telescopes)} telescopes"
            )

        alt_deg = alt_mean.to_value(u.deg)
        az_deg = az_mean.to_value(u.deg)

        if filename is None:
            filename = 'CTA-ULTRA6-LaPalma-div{}-az{}-alt{}.cfg'.format(
                str(self.div).replace(".", "_"),
                str(az_deg).replace(".", "_"),
                str(alt_deg).replace(".", "_"))

        path = Path(outdir) / filename

        with open(path, 'w') as f:
            f.write('#ifndef TELESCOPE\n')
            f.write('#  define TELESCOPE 0\n')
            f.write('#endif\n')
            f.write('#if TELESCOPE == 0\n')
            f.write('   TELESCOPE_THETA={:.2f} \n'.format(90 - alt_deg))
            f.write('   TELESCOPE_PHI={:.2f} \n'.format((360 - az_deg) % 360))
            f.write('\n% Global and default configuration for things missing in telescope-specific config.\n')
            f.write('#  include <{}>\n'.format(tel_configs[0]))
            for n, (tel, cfg) in enumerate(zip(self.telescopes, tel_configs), 1):
                f.write('\n#elif TELESCOPE == {:d}\n'.format(n))
                f.write('#  include <{}>\n'.format(cfg))
                f.write('   TELESCOPE_THETA={:.2f}\n'.format(
                    90 - tel.alt.to_value(u.deg)))
                # Az is unwrapped in divtel (arctan2 gives -180..180), so the
                # flip to sim_telarray's sense has to be brought back into
                # 0..360 rather than left negative.
                f.write('   TELESCOPE_PHI={:.2f}\n'.format(
                    (360 - tel.az.to_value(u.deg)) % 360))
            f.write('#else\n')
            f.write('   Error Invalid telescope for CTA-ULTRA6 La Palma configuration.\n')
            f.write('#endif\n')
            f.write('trigger_telescopes = 2 % In contrast to Prod-3 South we apply loose stereo trigger immediately\n')
            f.write('array_trigger = array_trigger_ultra6_diver-test.dat\n')

        if verbose:
            print(path.read_text())

        return path

    def display_2d(self, projection='xy', ax=None, **kwargs):
        """
        Display the array

        Parameters
        ----------
        projection: str
            'xy', 'xz' or 'yz'
        ax: `matplotlib.pyplot.axes`
        kwargs: args for `pyplot.scatter`

        Returns
        -------
        ax: `matplotlib.pyplot.axes`
        """
        ax = plt.gca() if ax is None else ax
        if 'color' not in kwargs:
            kwargs['color'] = 'black'

        # Axes are labelled in metres, so hand matplotlib bare numbers.
        positions_array = self.positions_array.to_value(u.m)
        barycenter = self.barycenter.to_value(u.m)

        if projection == 'xy':
            xx = positions_array[:, 1]
            yy = positions_array[:, 0]
            xb = barycenter[1]
            yb = barycenter[0]
            xv = self.pointing_vectors[:, 1]
            yv = self.pointing_vectors[:, 0]
            xlabel = 'y [m]'
            ylabel = 'x [m]'

        elif projection == 'xz':
            xx = positions_array[:, 0]
            yy = positions_array[:, 2]
            xb = barycenter[0]
            yb = barycenter[2]
            xv = self.pointing_vectors[:, 0]
            yv = self.pointing_vectors[:, 2]
            xlabel = 'x [m]'
            ylabel = 'z [m]'

        elif projection == 'yz':
            xx = positions_array[:, 1]
            yy = positions_array[:, 2]
            xb = barycenter[1]
            yb = barycenter[2]
            xv = self.pointing_vectors[:, 1]
            yv = self.pointing_vectors[:, 2]
            xlabel = 'y [m]'
            ylabel = 'z [m]'

        else:
            raise ValueError(f"projection should be either 'xy', ' yz' or 'xz' but is {projection}")

        scale = np.max(np.abs([xx, yy])) / 10.

        ax.scatter(xx, yy, **kwargs, label='telescopes')
        ax.scatter(xb, yb, marker='+', label='barycenter')
        ax.quiver(xx, yy, xv, yv,
                  color=kwargs['color'],
                  scale=scale,
                  )

        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.grid('on')
        # Pad through the margins rather than by setting limits outright:
        # `axis('equal')` keeps the axes box and stretches the data limits to
        # square up the aspect, so pinning both limits afterwards leaves it
        # nothing to stretch and it drops the y limits with a warning.
        ax.margins(0.25)
        ax.axis('equal')
        return ax

    def display_3d(self):
        """
        Display the array in 3d, with an arrow per telescope pointing

        Returns
        -------
        ax: `matplotlib.pyplot.axes`
        """
        # TODO: fix pointing quiver length issue
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        positions_array = self.positions_array.to_value(u.m)
        x = positions_array[:, 0]
        y = positions_array[:, 1]
        z = positions_array[:, 2]
        max_range = np.array([x.max() - x.min(), y.max() - y.min(), z.max() - z.min()]).max()
        scale = 1
        xb = scale * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + scale * (x.max() + x.min())
        yb = scale * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + scale * (y.max() + y.min())
        zb = scale * max_range * np.mgrid[-0.01:2:2, -0.01:2:2, -0.01:2:2][2].flatten() + scale * (z.max() + z.min())
        # Comment or uncomment following both lines to test the fake bounding box:
        for _xb, _yb, _zb in zip(xb, yb, zb):
            ax.plot([_xb], [_yb], [_zb], 'w')

        ax.quiver(x, y, z,
                  self.pointing_vectors[:, 0],
                  self.pointing_vectors[:, 1],
                  self.pointing_vectors[:, 2],
                  color='black',
                  length=max_range,
                  )

        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('z [m]')

        return ax

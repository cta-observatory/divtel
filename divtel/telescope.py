import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
    """

    _id = 0

    def __init__(self, x, y, z, focal, camera_radius):

        self.x = x.to(u.m)
        self.y = y.to(u.m)
        self.z = z.to(u.m)
        self.focal = focal.to(u.m)
        self.camera_radius = camera_radius.to(u.m)
        self.alt = u.Quantity(0, u.rad)
        self.az = u.Quantity(0, u.rad)
        Telescope._id += 1
        self.id = Telescope._id

    def point_to_altaz(self, alt, az):
        self.alt = alt.to(u.rad)
        self.az = az.to(u.rad)

    @property
    def zenith(self):
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
        return np.array([self.x.to(u.m).value, self.y.to(u.m).value, self.z.to(u.m).value]*u.m)

    def point_to_object(self, object):
        """
        Point to object.

        Parameters
        ----------
        object: numpy.array([x, y, z])
        """
        GT = np.sqrt(((self.position - object) ** 2).sum())
        alt_tel = np.arcsin((-self.z.value + object[2]) / GT)
        # Az runs clock-wise from X towards Y (see divtel.pointing), so the
        # Y offset enters negated: az = arctan2(-dy, dx).
        az_tel = np.arctan2(self.y.value - object[1], object[0] - self.x.value)
        self.point_to_altaz(alt_tel * u.rad, az_tel * u.rad)

    @property
    def pointing_vector(self):
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

    @property
    def positions_array(self):
        """
        All telescopes positions

        Returns
        -------
        positions_array: `numpy.array`
            2D numpy array
        """
        return np.array([tel.position for tel in self.telescopes])

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

        if div == 0:
            for tel in self.telescopes:
                tel.point_to_altaz(alt_mean, az_mean)
        else:
            g_point = pointing.pointG_position(self.barycenter, div, alt_mean, az_mean)
            for tel in self.telescopes:
                alt_tel, az_tel = pointing.tel_div_pointing(tel.position, g_point)
                tel.point_to_altaz(alt_tel*u.rad, az_tel*u.rad)

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

        if projection == 'xy':
            xx = self.positions_array[:, 1]
            yy = self.positions_array[:, 0]
            xb = self.barycenter[1]
            yb = self.barycenter[0]
            xv = self.pointing_vectors[:, 1]
            yv = self.pointing_vectors[:, 0]
            xlabel = 'y [m]'
            ylabel = 'x [m]'

        elif projection == 'xz':
            xx = self.positions_array[:, 0]
            yy = self.positions_array[:, 2]
            xb = self.barycenter[0]
            yb = self.barycenter[2]
            xv = self.pointing_vectors[:, 0]
            yv = self.pointing_vectors[:, 2]
            xlabel = 'x [m]'
            ylabel = 'z [m]'

        elif projection == 'yz':
            xx = self.positions_array[:, 1]
            yy = self.positions_array[:, 2]
            xb = self.barycenter[1]
            yb = self.barycenter[2]
            xv = self.pointing_vectors[:, 1]
            yv = self.pointing_vectors[:, 2]
            xlabel = 'y [m]'
            ylabel = 'z [m]'

        else:
            raise ValueError(f"projection should be either 'xy', ' yz' or 'xz' but is {projection}")

        scale = np.max([xx, yy]) / 10.

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
        # TODO: fix pointing quiver length issue
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        x = self.positions_array[:, 0]
        y = self.positions_array[:, 1]
        z = self.positions_array[:, 2]
        max_range = np.array([x.max() - x.min(), y.max() - y.min(), z.max() - z.min()]).max()
        scale = 1
        xb = scale * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + scale * (x.max() + x.min())
        yb = scale * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + scale * (y.max() + y.min())
        zb = scale * max_range * np.mgrid[-0.01:2:2, -0.01:2:2, -0.01:2:2][2].flatten() + scale * (z.max() + z.min())
        # Comment or uncomment following both lines to test the fake bounding box:
        for xb, yb, zb in zip(xb, yb, zb):
            ax.plot([xb], [yb], [zb], 'w')

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

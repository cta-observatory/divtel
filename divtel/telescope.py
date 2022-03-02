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
        return u.Quantity(np.pi * (self.camera_radius / self.focal)**2, u.rad**2)

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
        az_tel = np.arctan2((- self.y.value + object[1]), (- self.x.value + object[0]))
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
        """
        if div == 0:
            for tel in self.telescopes:
                tel.point_to_altaz(alt_mean, az_mean)
        else:
            g_point = pointing.pointG_position(self.barycenter, div, alt_mean, az_mean)
            for tel in self.telescopes:
                alt_tel, az_tel = pointing.tel_div_pointing(tel.position, g_point)
                tel.point_to_altaz(alt_tel*u.rad, az_tel*u.rad)

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
        ax.axis('equal')
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.set_xlim(xlim[0] - 0.25 * np.abs(xlim[0]), xlim[1] + 0.25 * np.abs(xlim[1]))
        ax.set_ylim(ylim[0] - 0.25 * np.abs(ylim[0]), ylim[1] + 0.25 * np.abs(ylim[1]))
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

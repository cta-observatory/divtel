import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from mpl_toolkits.axisartist import Subplot
from mpl_toolkits.axisartist import SubplotHost, ParasiteAxesAuxTrans
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
import astropy.units as u


def display_array_pointing_in_sky(array):
    #TODO: take an array class and plot the pointings FoV in the sky
    # need to handle projection
    raise NotImplementedError("TODO")


def sky_fov(telescope, ax=None):
    """
    Display the telescope FoV in the sky

    Parameters
    ----------
    telescope: `Telescope`
    ax: `matplotlib.pyplot.axes`

    Returns
    -------
    ax: `matplotlib.pyplot.axes`
    """
    raise NotImplementedError("TODO")


def polar_stuff(fig, telescope):
    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
    tr = Affine2D().scale(np.pi/180., 1.).translate(+np.pi/2.,0) + PolarAxes.PolarTransform()

    # polar projection, which involves cycle, and also has limits in
    # its coordinates, needs a special method to find the extremes
    # (min, max of the coordinate within the view).

    # 20, 20 : number of sampling points along x, y direction
    n = 1
    extreme_finder = angle_helper.ExtremeFinderCycle(n, n,
                                                     lon_cycle=360,
                                                     lat_cycle=None,
                                                     lon_minmax=None,
                                                     lat_minmax=(-90, 90),
                                                     )

    grid_locator1 = angle_helper.LocatorDMS(12)
    # Find a grid values appropriate for the coordinate (degree,
    # minute, second).

    tick_formatter1 = angle_helper.FormatterDMS()
    # And also uses an appropriate formatter.  Note that,the
    # acceptable Locator and Formatter class is a bit different than
    # that of mpl's, and you cannot directly use mpl's Locator and
    # Formatter here (but may be possible in the future).

    grid_helper = GridHelperCurveLinear(tr,
                                        extreme_finder=extreme_finder,
                                        grid_locator1=grid_locator1,
                                        tick_formatter1=tick_formatter1
                                        )

    ax1 = SubplotHost(fig, 1, 1, 1, grid_helper=grid_helper)

    # make ticklabels of right and top axis visible.
    ax1.axis["right"].major_ticklabels.set_visible(True)
    ax1.axis["top"].major_ticklabels.set_visible(True)

    # let right axis shows ticklabels for 1st coordinate (angle)
    ax1.axis["right"].get_helper().nth_coord_ticks = 0
    # let bottom axis shows ticklabels for 2nd coordinate (radius)
    ax1.axis["bottom"].get_helper().nth_coord_ticks = 1

    fig.add_subplot(ax1)

    # A parasite axes with given transform
    ax2 = ParasiteAxesAuxTrans(ax1, tr, "equal")
    # note that ax2.transData == tr + ax1.transData
    # Anything you draw in ax2 will match the ticks and grids of ax1.
    ax1.parasites.append(ax2)
    # intp = cbook.simple_linear_interpolation
    #ax2.plot(intp(np.array([0, 30]), 50),
    #         intp(np.array([10., 10.]), 50),
    #         linewidth=2.0)



    x = np.rad2deg(telescope.az.value) * np.cos(telescope.alt.value)
    y = np.rad2deg(telescope.alt.value)

    circle = plt.Circle((np.rad2deg(telescope.az.value - np.pi) * np.sin(telescope.alt.value),
                         np.rad2deg(-telescope.alt.value * np.cos((telescope.az.value - np.pi)))),
                        radius=7.7 / 2,
                        color="red",
                        alpha=0.2,
                        )

    circle = plt.Circle((x, y),
                        radius=7.7 / 2,
                        color="red",
                        alpha=0.2,
                        )
    ax1.add_artist(circle)
    # point = ax1.scatter(x, y, c="b", s=20, zorder=10, transform=ax2.transData)
    ax2.annotate(1, (x, y), fontsize=15, xytext=(4, 4), textcoords='offset pixels')

    ax1.set_xlim(-180, 180)
    ax1.set_ylim(0, 90)
    ax1.set_aspect(1.)
    ax1.grid(True, zorder=0)
    ax1.set_xlabel("Azimuth in degrees", fontsize=20)
    ax1.set_ylabel("Zenith in degrees", fontsize=20)

    plt.show()
    return fig




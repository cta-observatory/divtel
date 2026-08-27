import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from mpl_toolkits.axisartist import HostAxes
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
import astropy.units as u


def display_hyper_fov(array, ax=None, m_cut=1, cmap="viridis", show_area=True):
    """
    Display the array's hyper field of view on the sky.

    Each telescope sees a disc; divergent pointing spreads those discs out,
    widening the sky the array covers at the cost of how many telescopes see
    any given part of it. This draws the union of the discs, shaded by
    multiplicity -- the number of telescopes seeing each patch. Patches seen by
    two or more can be reconstructed stereoscopically; patches seen by one
    cannot.

    Parameters
    ----------
    array: `Array`
    ax: `matplotlib.pyplot.axes`, optional
    m_cut: int
        multiplicity counted towards the area reported in the title; patches
        below it are still drawn, faded.
    cmap: str
        colormap used to shade multiplicity
    show_area: bool
        title the axes with the covered area

    Returns
    -------
    ax: `matplotlib.pyplot.axes`

    Notes
    -----
    Drawn in an equal-area projection centred on the array's mean pointing, so
    the axes are degrees of offset from that direction and the areas are true
    solid angles -- see `Array.hyper_fov`.
    """
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import BoundaryNorm
    from matplotlib.patches import Polygon as MplPolygon

    ax = plt.gca() if ax is None else ax

    area, patches = array.hyper_fov(m_cut=m_cut)
    if not patches:
        raise ValueError("the array covers no sky; are the telescopes pointed?")

    multiplicities = [m for _, m in patches]
    m_max = max(multiplicities)

    # One colour per integer multiplicity, so the shading is readable as a
    # count rather than a gradient.
    bounds = np.arange(0.5, m_max + 1.5)
    norm = BoundaryNorm(bounds, plt.get_cmap(cmap).N)

    polygons, colours = [], []
    for patch, multiplicity in patches:
        x, y = patch.exterior.xy
        polygons.append(MplPolygon(np.column_stack([x, y]), closed=True))
        colours.append(multiplicity)

    collection = PatchCollection(polygons, cmap=cmap, norm=norm,
                                 edgecolor="white", linewidth=0.4)
    collection.set_array(np.array(colours))
    # Patches below the cut are drawn but faded, so the shape of the whole
    # covered area stays visible while the counted part stands out.
    collection.set_alpha(None)
    ax.add_collection(collection)

    for patch, multiplicity in patches:
        if multiplicity < m_cut:
            x, y = patch.exterior.xy
            ax.fill(x, y, facecolor="white", alpha=0.55, edgecolor="none", zorder=2)

    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("offset from array pointing [deg]")
    ax.set_ylabel("offset from array pointing [deg]")
    ax.margins(0.08)

    colourbar = ax.figure.colorbar(collection, ax=ax,
                                   ticks=np.arange(1, m_max + 1),
                                   label="telescopes seeing this patch")
    colourbar.ax.set_yticklabels([str(i) for i in range(1, m_max + 1)])

    if show_area:
        label = "covered" if m_cut <= 1 else f"seen by $\\geq${m_cut}"
        ax.set_title(f"hyper FoV: {area.value:.1f} deg$^2$ {label}")

    return ax


def display_array_pointing_in_sky(array, ax=None, **kwargs):
    """
    Deprecated alias for `display_hyper_fov`.
    """
    return display_hyper_fov(array, ax=ax, **kwargs)


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

    # number of sampling points along the x and y directions; the grid finder
    # needs more than one or it sees a zero-width interval and divides by it.
    n = 20
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

    ax1 = fig.add_subplot(1, 1, 1, axes_class=HostAxes, grid_helper=grid_helper)

    # make ticklabels of right and top axis visible.
    ax1.axis["right"].major_ticklabels.set_visible(True)
    ax1.axis["top"].major_ticklabels.set_visible(True)

    # let right axis shows ticklabels for 1st coordinate (angle)
    ax1.axis["right"].get_helper().nth_coord_ticks = 0
    # let bottom axis shows ticklabels for 2nd coordinate (radius)
    ax1.axis["bottom"].get_helper().nth_coord_ticks = 1

    # A parasite axes with the given transform: ax2.transData == tr +
    # ax1.transData, so anything drawn in ax2 matches the ticks and grids of
    # ax1. get_aux_axes builds it and registers it as a parasite for us.
    ax2 = ax1.get_aux_axes(tr, viewlim_mode="equal")
    # intp = cbook.simple_linear_interpolation
    #ax2.plot(intp(np.array([0, 30]), 50),
    #         intp(np.array([10., 10.]), 50),
    #         linewidth=2.0)



    x = np.rad2deg(telescope.az.value) * np.cos(telescope.alt.value)
    y = np.rad2deg(telescope.alt.value)

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
    ax1.set_ylabel("Altitude in degrees", fontsize=20)

    plt.show()
    return fig




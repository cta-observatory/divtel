import matplotlib.pyplot as plt
import numpy as np


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

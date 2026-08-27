import astropy.units as u
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


def multiplicity_plot(array, m_cut=1, ax=None, cmap="viridis"):
    """
    Bar chart of how much sky is seen by how many telescopes.

    `display_hyper_fov` shades the sky map by multiplicity, which shows where
    the well-covered parts are; this counts them up, which shows how much there
    is of each. Bars are coloured to match that map, so the two read together.

    Parameters
    ----------
    array: `Array`
    m_cut: int
        multiplicity counted towards the area reported in the title; bars below
        it are still drawn, faded
    ax: `matplotlib.pyplot.axes`, optional
    cmap: str
        colormap used to shade multiplicity, as in `display_hyper_fov`

    Returns
    -------
    ax: `matplotlib.pyplot.axes`

    Examples
    --------
    >>> array.divergent_pointing(0.05, 70 * u.deg, 180 * u.deg)
    >>> multiplicity_plot(array)
    """
    from matplotlib.colors import BoundaryNorm

    area, patches = array.hyper_fov(m_cut=m_cut)
    if not patches:
        raise ValueError("the array covers no sky; are the telescopes pointed?")

    multiplicity, per_multiplicity = array.multiplicity_profile(patches=patches)
    mean, variance = array.multiplicity_moments(patches=patches)

    ax = plt.gca() if ax is None else ax

    # Same discrete colour scale as the sky map, so a bar and the patches it
    # counts come out the same colour. That means normalising over the
    # multiplicities actually present, exactly as `display_hyper_fov` does --
    # normalising over the telescope count instead would leave every bar of a
    # strongly divergent array at the dark end of the map.
    colors = plt.get_cmap(cmap)
    norm = BoundaryNorm(np.arange(0.5, multiplicity.max() + 1.5), colors.N)
    bars = ax.bar(multiplicity, per_multiplicity.to_value(u.deg ** 2),
                  color=[colors(norm(m)) for m in multiplicity],
                  edgecolor="black", linewidth=0.5)

    for bar, m in zip(bars, multiplicity, strict=True):
        if m < m_cut:
            bar.set_alpha(0.25)

    label = "covered" if m_cut <= 1 else rf"seen by $\geq${m_cut}"
    ax.set_title(f"hyper FoV: {area.value:.1f} deg$^2$ {label}\n"
                 rf"multiplicity {mean:.1f} $\pm$ {np.sqrt(variance):.1f}")
    ax.set_xlabel("multiplicity")
    ax.set_ylabel("hyper FoV [deg$^2$]")
    ax.set_xticks(multiplicity)
    return ax


def display_groups(groups, projection="xy", ax=None, cmap="tab10", annotate=True):
    """
    Display sub-arrays, each in its own colour with its own barycenter.

    `Array.display_2d` draws the whole array in one colour around one
    barycenter. On a layout of several telescope types that hides the thing
    worth seeing: where each type sits, and how far each type's mean pointing
    has swung away from the others under divergence.

    Parameters
    ----------
    groups: dict of str to `Array`
        as returned by `Array.group_by`
    projection: str
        'xy', 'xz' or 'yz'
    ax: `matplotlib.pyplot.axes`, optional
    cmap: str
        colormap the group colours are taken from
    annotate: bool
        label each barycenter with its group name

    Returns
    -------
    ax: `matplotlib.pyplot.axes`

    Examples
    --------
    >>> groups = array.group_by({"LST": range(1, 5), "MST": range(5, 20)})
    >>> display_groups(groups)
    """
    axes = {"xy": (1, 0), "xz": (0, 2), "yz": (1, 2)}
    if projection not in axes:
        raise ValueError(
            f"projection should be either 'xy', 'yz' or 'xz' but is {projection}"
        )
    if not groups:
        raise ValueError("no groups to display")

    horizontal, vertical = axes[projection]
    labels = "xyz"

    ax = plt.gca() if ax is None else ax
    colors = plt.get_cmap(cmap)

    # One scale for every group, so arrow lengths stay comparable between them.
    # `scale_units="xy"` puts arrow length in data units, which is the only way
    # to size them against the array: a unit pointing vector is drawn a tenth of
    # the array across. Note that an arrow shows the pointing *projected into
    # this plane*, so it shortens as the array points out of it.
    everything = np.concatenate(
        [array.positions_array.to_value(u.m) for array in groups.values()]
    )
    extent = max(np.ptp(everything[:, horizontal]), np.ptp(everything[:, vertical]))
    scale = 10.0 / extent if extent > 0 else 1.0
    quiver_style = {"angles": "xy", "scale_units": "xy", "scale": scale}

    for i, (name, array) in enumerate(groups.items()):
        color = colors(i % colors.N)
        positions = array.positions_array.to_value(u.m)
        barycenter = array.barycenter.to_value(u.m)
        pointings = array.pointing_vectors

        xx, yy = positions[:, horizontal], positions[:, vertical]
        ax.scatter(xx, yy, color=color, label=name)
        ax.quiver(xx, yy, pointings[:, horizontal], pointings[:, vertical],
                  color=color, alpha=0.4, **quiver_style)

        xb, yb = barycenter[horizontal], barycenter[vertical]
        ax.scatter(xb, yb, marker="+", s=200, linewidths=2, color=color)
        ax.quiver(xb, yb,
                  pointings[:, horizontal].mean(), pointings[:, vertical].mean(),
                  color=color, **quiver_style)

        if annotate:
            ax.annotate(name, (xb, yb), textcoords="offset points",
                        xytext=(8, 8), color=color, fontweight="bold")

    ax.set_xlabel(f"{labels[horizontal]} [m]")
    ax.set_ylabel(f"{labels[vertical]} [m]")
    ax.grid("on")
    ax.legend(frameon=False)
    ax.margins(0.25)
    ax.axis("equal")
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

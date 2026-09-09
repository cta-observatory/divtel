import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np


def display_hyper_fov(array, ax=None, min_telescopes=2, cmap="viridis", show_area=True):
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
    min_telescopes: int
        multiplicity counted towards the area reported in the title; patches
        below it are still drawn, faded. The default, 2, reports the
        stereoscopic area; pass 1 for the whole covered area.
    cmap: str
        colormap used to shade multiplicity
    show_area: bool
        title the axes with the covered area

    Returns
    -------
    ax: `matplotlib.pyplot.axes`

    Notes
    -----
    Drawn in an equal-area projection centred on the array's mean pointing:
    x is offset in azimuth, y is offset in altitude, and the areas are true
    solid angles -- see `Array.hyper_fov`.
    """
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import BoundaryNorm
    from matplotlib.patches import Polygon as MplPolygon

    ax = plt.gca() if ax is None else ax

    area, patches = array.hyper_fov(min_telescopes=min_telescopes)
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
        if multiplicity < min_telescopes:
            x, y = patch.exterior.xy
            ax.fill(x, y, facecolor="white", alpha=0.55, edgecolor="none", zorder=2)

    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("azimuth offset [deg]")
    ax.set_ylabel("altitude offset [deg]")
    ax.margins(0.08)

    colourbar = ax.figure.colorbar(collection, ax=ax,
                                   ticks=np.arange(1, m_max + 1),
                                   label="telescopes seeing this patch")
    colourbar.ax.set_yticklabels([str(i) for i in range(1, m_max + 1)])

    if show_area:
        label = "covered" if min_telescopes <= 1 else f"seen by $\\geq${min_telescopes}"
        ax.set_title(f"hyper FoV: {area.value:.1f} deg$^2$ {label}")

    return ax


def multiplicity_plot(array, min_telescopes=2, ax=None, cmap="viridis"):
    """
    Bar chart of how much sky is seen by how many telescopes.

    `display_hyper_fov` shades the sky map by multiplicity, which shows where
    the well-covered parts are; this counts them up, which shows how much there
    is of each. Bars are coloured to match that map, so the two read together.

    Parameters
    ----------
    array: `Array`
    min_telescopes: int
        multiplicity counted towards the area reported in the title; bars below
        it are still drawn, faded. The default, 2, reports the stereoscopic
        area; pass 1 for the whole covered area.
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

    area, patches = array.hyper_fov(min_telescopes=min_telescopes)
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
        if m < min_telescopes:
            bar.set_alpha(0.25)

    label = "covered" if min_telescopes <= 1 else rf"seen by $\geq${min_telescopes}"
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


class Projection:
    """
    A Lambert azimuthal equal-area frame, centred on a region.

    Areas on the page are true solid angles -- the same projection
    `Array.hyper_fov` measures in, chosen for the same reasons: it keeps areas
    true and it has no pole to fall over near the zenith.

    Held as an object rather than recomputed per call so that everything drawn
    on one panel shares it. That matters as soon as more than the region itself
    is drawn: credible bands cut at different levels have different centroids,
    and projecting each in its own frame would slide them apart instead of
    nesting them, while camera footprints in a frame of their own would land
    somewhere other than where they point.

    Build one with `projection_frame` and pass it to `project_region`,
    `project_directions` and `camera_rims`.
    """

    def __init__(self, centre, east, north, offset=(0.0, 0.0), turn=None):
        self.centre = centre
        self.east = east
        self.north = north
        self.offset = offset
        self.turn = turn

    def __call__(self, directions):
        """
        Project unit vectors onto the page.

        Parameters
        ----------
        directions: `numpy.ndarray`
            shape (n, 3), unit vectors in the region's frame

        Returns
        -------
        x, y: `numpy.ndarray`
            degrees from the frame's centre
        """
        directions = np.atleast_2d(directions)
        along = np.clip(directions @ self.centre, -1.0, 1.0)
        scale = np.degrees(np.sqrt(2.0 / np.maximum(1.0 + along, 1e-12)))
        x = scale * (directions @ self.east) - self.offset[0]
        y = scale * (directions @ self.north) - self.offset[1]
        if self.turn is None:
            return x, y
        return tuple((np.column_stack([x, y]) @ self.turn.T).T)


def projection_frame(region, align=True):
    """
    The equal-area frame a region is drawn in.

    Parameters
    ----------
    region: `divtel.region.SkyRegion`
        placed over the array. The frame is centred on its centroid.
    align: bool
        turn the region so its long axis lies along x. The rotation is rigid
        and happens in the equal-area plane, so areas and camera shapes are
        untouched; only the frame the reader sees turns, which is why the axes
        stop being azimuth and altitude when this is on.

    Returns
    -------
    `Projection`
    """
    from .pointing import as_altaz, local_frame

    centre = region.centroid
    east, north = local_frame(*as_altaz(centre))
    frame = Projection(centre, east, north)

    if not align:
        return frame

    weights = region.weights
    x, y = frame(region.directions)
    offset = (float(np.average(x, weights=weights)),
              float(np.average(y, weights=weights)))

    points = np.column_stack([x - offset[0], y - offset[1]])
    spread = (points * weights[:, None]).T @ points
    principal = np.linalg.eigh(spread)[1][:, -1]
    turn = np.array([[principal[0], principal[1]],
                     [-principal[1], principal[0]]])
    return Projection(centre, east, north, offset, turn)


def project_region(region, align=True, frame=None):
    """
    A placed region as x and y on the page, in degrees.

    Shared rather than inlined because a figure's *size* has to be chosen from
    the same projection its *contents* are drawn in, and the two drifting apart
    is not a visible failure -- it is a panel that is quietly the wrong shape.

    Parameters
    ----------
    region: `divtel.region.SkyRegion`
        placed over the array
    align: bool
        see `projection_frame`; ignored when `frame` is given
    frame: `Projection`, optional
        draw in this frame instead of the region's own. Pass one region's frame
        when projecting another -- a narrower credible band, say -- so the two
        land on top of each other.

    Returns
    -------
    x, y: `numpy.ndarray`
        degrees from the frame's centre
    """
    frame = frame if frame is not None else projection_frame(region, align)
    return frame(region.directions)


def project_directions(directions, frame):
    """
    Arbitrary unit vectors on the page, in a region's frame.

    Parameters
    ----------
    directions: `numpy.ndarray`
        shape (n, 3)
    frame: `Projection`

    Returns
    -------
    x, y: `numpy.ndarray`
    """
    return frame(directions)


def region_span(region, align=True):
    """
    Height over width of a region on the page, for choosing a figure's shape.

    Panels drawn by `multiplicity_over_region` are equal-aspect, so a figure not
    shaped like the region leaves its axes shrunk inside their slots with the
    colour bars standing beside empty space.

    Parameters
    ----------
    region: `divtel.region.SkyRegion`
    align: bool
        must match what the panels are drawn with

    Returns
    -------
    float
    """
    x, y = project_region(region, align=align)
    return float(np.ptp(y) / max(np.ptp(x), 1e-9))


def type_colors(array, palette=("#2a78d6", "#eb6834", "#1baf7a", "#eda100"),
                names=None):
    """
    A colour per telescope, by the angular radius of its camera.

    What every pointing strategy means by a telescope's type -- the same split
    `Array.group_by("fov_radius")` makes -- so a figure coloured this way shows
    which instrument is being sent where. Free pointing puts the narrow cameras
    where the map is sharp and the wide ones on the tails, which is a degree of
    freedom neither divergence nor an even sub-array split can use, and it is
    invisible without the colours.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    palette: sequence of str
        colours to draw from, widest camera first
    names: sequence of str, optional
        one name per camera type, widest first, used in the legend instead of
        the angular radius (e.g. ``("SST", "MST")`` for CTAO-South). Left out,
        the legend names each group by its field-of-view radius, which is all
        `fov_radius` knows -- it has no notion of "SST" or "MST".

    Returns
    -------
    colors: list of str
        one per telescope, in ``array.telescopes`` order
    labels: dict
        colour to a legend label for that camera type
    """
    radii = sorted({t.fov_radius for t in array.telescopes}, reverse=True)
    if names is not None and len(names) != len(radii):
        raise ValueError(
            f"names has {len(names)} entries but the array has {len(radii)} "
            "camera types")
    by_radius = {radius: palette[index % len(palette)]
                 for index, radius in enumerate(radii)}
    if names is None:
        labels = {colour: f"{radius.to_value(u.deg):.1f}$\\degree$ camera"
                  for radius, colour in by_radius.items()}
    else:
        labels = {by_radius[radius]: f"{name} camera"
                  for radius, name in zip(radii, names)}
    return [by_radius[t.fov_radius] for t in array.telescopes], labels


def camera_rims(array, region=None, ax=None, frame=None, colors=None,
                labels=None, type_names=None, points=129, linewidth=1.0,
                alpha=0.75, zorder=6):
    """
    The outline of every telescope's field of view, drawn over a region.

    The picture behind a coverage number: which parts of the sky the cameras
    actually sit on, and which parts they miss. Composes with
    `multiplicity_over_region` and `probability_over_region` -- pass the same
    `frame` to all of them and they land on top of each other.

    Each rim is sampled on the sphere and then projected, rather than drawn as
    a circle on the page. Away from the centre of an equal-area projection the
    image of a circle is not one, and it is the true shape that shows a camera
    at the edge of the frame covering what it really covers.

    Parameters
    ----------
    array: `divtel.telescope.Array`
        already pointed, by any strategy
    region: `divtel.region.SkyRegion`, optional
        used only to build a frame when none is given
    ax: `matplotlib.axes.Axes`, optional
    frame: `Projection`, optional
        the frame to draw in; built from `region` otherwise
    colors: list of str or ``"type"``, optional
        one colour per telescope, in ``array.telescopes`` order. ``"type"``
        colours by camera radius via `type_colors`, which is usually what is
        wanted. A single neutral ink otherwise.
    labels: dict, optional
        colour to legend label. Filled in by ``colors="type"``.
    type_names: sequence of str, optional
        passed to `type_colors` when ``colors="type"``, to legend the camera
        types by name (e.g. ``("SST", "MST")``) rather than by radius.
    points: int
        samples around each rim
    linewidth, alpha, zorder: float
        passed to the line

    Returns
    -------
    extent: tuple of float
        ``(x_min, x_max, y_min, y_max)`` of the rims on the page, in degrees.
        Frame the panel on this together with the region's own extent, or a
        camera hanging off the edge is silently cropped and the picture
        overstates the coverage. `frame_on` does that.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6.4, 5.2))
    if frame is None:
        if region is None:
            raise ValueError("camera_rims needs a region or a frame to draw in")
        frame = projection_frame(region)

    if colors == "type":
        colors, labels = type_colors(array, names=type_names)
    elif colors is None:
        colors = ["#52514e"] * len(array.telescopes)

    angles = np.linspace(0, 2 * np.pi, points)
    extent = None
    for telescope, direction, colour in zip(array.telescopes,
                                            array.pointing_vectors, colors,
                                            strict=True):
        radius = telescope.fov_radius.to_value(u.rad)
        # Two directions perpendicular to where this telescope points, so the
        # rim can be swept around it. Any pair will do; taking one from the
        # frame's centre keeps the sampling steady as the pointing moves.
        first = np.cross(direction, frame.centre)
        norm = np.linalg.norm(first)
        first = frame.east.copy() if norm < 1e-9 else first / norm
        second = np.cross(direction, first)

        rim = (np.cos(radius) * direction
               + np.sin(radius) * (np.cos(angles)[:, None] * first
                                   + np.sin(angles)[:, None] * second))
        x, y = frame(rim)
        extent = _union(extent, (x.min(), x.max(), y.min(), y.max()))
        ax.plot(x, y, linewidth=linewidth, color=colour, alpha=alpha,
                zorder=zorder)

    for colour, label in (labels or {}).items():
        ax.plot([], [], linewidth=1.6, color=colour, label=label)

    return extent


def _union(*extents):
    """The smallest (x_min, x_max, y_min, y_max) containing all of them."""
    boxes = [box for box in extents if box is not None]
    return (min(b[0] for b in boxes), max(b[1] for b in boxes),
            min(b[2] for b in boxes), max(b[3] for b in boxes))


def region_extent(region, frame=None, align=True):
    """
    Where a region falls on the page.

    Parameters
    ----------
    region: `divtel.region.SkyRegion`
    frame: `Projection`, optional
    align: bool

    Returns
    -------
    tuple of float
        ``(x_min, x_max, y_min, y_max)``, in degrees
    """
    x, y = project_region(region, align=align, frame=frame)
    return float(x.min()), float(x.max()), float(y.min()), float(y.max())


def frame_on(ax, *extents, pad=0.04, equal=True):
    """
    Set a panel's limits to hold everything drawn on it.

    A region and the cameras pointed at it need framing together: on the
    region alone a camera at the edge is silently cropped and the picture
    overstates the coverage, and on the cameras alone a strategy that flings
    them wide shrinks the region to a smear.

    Parameters
    ----------
    ax: `matplotlib.axes.Axes`
    *extents: tuple of float
        ``(x_min, x_max, y_min, y_max)`` boxes, as `region_extent` and
        `camera_rims` return them
    pad: float
        margin to leave, as a fraction of the box
    equal: bool
        keep the aspect equal, so areas on the page stay true areas on the sky

    Returns
    -------
    `matplotlib.axes.Axes`
    """
    x_min, x_max, y_min, y_max = _union(*extents)
    margin = pad * max(x_max - x_min, y_max - y_min)
    ax.set_xlim(x_min - margin, x_max + margin)
    ax.set_ylim(y_min - margin, y_max + margin)
    if equal:
        ax.set_aspect("equal")
    return ax


def probability_over_region(region, bands=None, ax=None, frame=None,
                            color="#2a78d6", size=4.0, align=True):
    """
    A region shaded by where the source probably is.

    The companion to `multiplicity_over_region`: that one shows what the array
    delivers, this one shows what it is being asked to cover. Draw the two side
    by side, or draw `camera_rims` over this one to see which parts of the
    probability the cameras sit on.

    Shaded by credible level, never by raw probability density. Density spans
    orders of magnitude and plots as a thin bright thread with nothing around
    it; the credible level answers the question actually being asked, which is
    where the source is likely to be and how confident that is.

    With no `bands` the level is computed from the region's own weights, the
    same greedy construction that cut it: take the most probable direction,
    then the next, and shade each by the probability accumulated when it was
    taken. That gives a continuous ramp rather than a few discrete rings, and
    it needs nothing but the region.

    Parameters
    ----------
    region: `divtel.region.SkyRegion`
        placed over the array. Sets the frame, and is what gets shaded when
        `bands` is None.
    bands: dict of float to `divtel.region.SkyRegion`, optional
        nested credible regions, keyed by level, drawn as discrete rings
        instead. **Place them with one shared rotation** --
        `divtel.region.SkyRegion.placement_rotation` and `rotate` -- or they
        will not nest.
    ax: `matplotlib.axes.Axes`, optional
    frame: `Projection`, optional
    color: str
        the hue the bands are shaded in, light to dark
    size: float
        marker area
    align: bool
        see `projection_frame`; ignored when `frame` is given

    Returns
    -------
    frame: `Projection`
        the frame everything was drawn in, to pass to `camera_rims`
    """
    from matplotlib.colors import LinearSegmentedColormap, to_rgb

    if ax is None:
        _, ax = plt.subplots(figsize=(6.4, 5.2))
    if frame is None:
        frame = projection_frame(region, align)

    if bands is None:
        # The greedy credible-level construction, on the region's own weights.
        weights = region.weights / region.weights.sum()
        order = np.argsort(weights)[::-1]
        level = np.empty(len(weights))
        level[order] = np.cumsum(weights[order])

        x, y = frame(region.directions)
        base = to_rgb(color)
        ramp = LinearSegmentedColormap.from_list(
            "credible", [tuple(1 - (1 - c) * w for c in base)
                         for w in np.linspace(1.0, 0.18, 32)])
        ax.scatter(x, y, c=level, s=size, marker=".", linewidths=0, cmap=ramp,
                   vmin=0.0, vmax=float(level.max()), zorder=2)

        ax.set_xlabel("along the region [deg]" if frame.turn is not None
                      else "offset in azimuth [deg]")
        ax.set_ylabel("across it [deg]" if frame.turn is not None
                      else "offset in altitude [deg]")
        ax.set_aspect("equal")
        return frame

    levels = sorted(bands)

    base = to_rgb(color)
    # One hue, light to dark: the bands are a magnitude, not several identities.
    shades = dict(zip(levels,
                      [tuple(1 - (1 - c) * w for c in base)
                       for w in np.linspace(0.30, 1.0, len(levels))[::-1]],
                      strict=True))

    # Outermost first, so the darker inner bands land on top rather than under.
    # The legend entries are added in the other order, so a reader meets them
    # innermost first, the way the bands are usually quoted.
    for depth, level in enumerate(reversed(levels)):
        x, y = frame(bands[level].directions)
        ax.scatter(x, y, s=size, marker=".", linewidths=0, color=shades[level],
                   zorder=2 + depth)
    for level in levels:
        ax.scatter([], [], s=26, marker=".", linewidths=0, color=shades[level],
                   label=f"{level:.0%} credible")

    ax.set_xlabel("along the region [deg]" if frame.turn is not None
                  else "offset in azimuth [deg]")
    ax.set_ylabel("across it [deg]" if frame.turn is not None
                  else "offset in altitude [deg]")
    ax.set_aspect("equal")
    return frame


def multiplicity_over_region(array, region, ax=None, title=None, vmax=None,
                             radius=None, colorbar=True, align=True, size=4.0,
                             cmap="Blues", blind_color="0.85", frame=None):
    """
    A weighted region coloured by how many telescopes see each part of it.

    The figure that separates the pointing strategies. A parallel array paints
    one solid blob and leaves the rest blank; a diverged one paints a wide ring,
    evenly; a shaped one paints the region's own shape, darkest where the weight
    is. Whether it does is something to look at rather than take on trust from a
    single number.

    Both the region and the colours are in the equal-area projection of
    `project_region`, so a wide patch of low multiplicity looks as big as it is.

    Parameters
    ----------
    array: `divtel.telescope.Array`
        already pointed, by any strategy
    region: `divtel.region.SkyRegion`
        placed over the array
    ax: `matplotlib.axes.Axes`, optional
    title: str, optional
    vmax: int, optional
        top of the colour scale, above which the colour stops changing. Set it
        by hand across a row of panels, or each is normalised to its own maximum
        and they cannot be compared. A parallel array piles every telescope on
        one spot and its maximum is the whole array, so scaling a row to *that*
        leaves every other panel a uniform pale wash: pick a high percentile of
        the row instead and let the parallel panel clip.
    radius: `astropy.Quantity`, optional
        half-width of the plot; sized to hold the region by default
    colorbar: bool
        draw one beside this panel. Turn it off for a row of panels sharing a
        scale and add a single bar for the figure.
    align: bool
        turn the region so its long axis lies along the page. A region that runs
        diagonally fills a square panel with empty corners and packs badly beside
        other panels; laid flat the same region is wide and short. The axis
        labels stop being azimuth and altitude when this is on, which is why it
        can be turned off -- for a compact region there is nothing to align and
        the true orientation is worth more. Ignored when `frame` is given.
    size: float
        marker area for the directions
    cmap: str
        one hue, light to dark: multiplicity is a magnitude, not an identity
    blind_color: str
        colour for directions no telescope sees
    frame: `Projection`, optional
        draw in this frame rather than the region's own, so this panel lands on
        top of whatever else was drawn in it -- `camera_rims`, say

    Returns
    -------
    `matplotlib.collections.PathCollection`
        the drawn points, to hang a shared colour bar on

    Notes
    -----
    Directions no telescope sees are drawn in grey rather than left off the page.
    Where an array is blind is most of the point of the picture, and a gap reads
    as "no region here" instead of "no telescopes here".
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6.4, 5.2))

    if frame is None:
        frame = projection_frame(region, align)
    multiplicity = region.multiplicity(array)
    x, y = frame(region.directions)

    turned = frame.turn is not None
    ax.set_xlabel("along the region [deg]" if turned else "offset in azimuth [deg]")
    ax.set_ylabel("across it [deg]" if turned else "offset in altitude [deg]")
    if title:
        ax.set_title(title, loc="left")

    blind = multiplicity < 1
    if blind.any():
        ax.scatter(x[blind], y[blind], s=size, marker=".", linewidths=0,
                   color=blind_color, zorder=1)

    seen = ~blind
    # The scale starts at 1, not 0: an unseen direction is grey and out of the
    # ramp altogether, so giving zero a colour would waste the pale end of it on
    # a value that never appears.
    dots = ax.scatter(x[seen], y[seen], c=multiplicity[seen], s=size, marker=".",
                      linewidths=0, cmap=cmap, vmin=1,
                      vmax=vmax if vmax is not None else max(2, multiplicity.max()),
                      zorder=2)

    if colorbar:
        bar = ax.figure.colorbar(dots, ax=ax, fraction=0.045, pad=0.02)
        bar.set_label("telescopes seeing this direction")

    if radius is not None:
        half = radius.to_value(u.deg)
        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)
    ax.set_aspect("equal")
    return dots


def multiplicity_by_probability(profiles, ax=None, m_cut=2, bins=5):
    """
    Mean multiplicity in each fifth of the weight, least likely first.

    The plain reading of whether an array is looking hardest where the source
    probably is: a flat set of bars is an array spread evenly over the region, a
    rising one is `divtel.strategy.shaped_pointing` doing what it claims.

    Parameters
    ----------
    profiles: dict of str to list of float
        one entry per strategy, as `divtel.strategy.multiplicity_by_probability`
        returns them
    ax: `matplotlib.axes.Axes`, optional
    m_cut: int
        drawn as a floor line -- below it a direction is covered but not usable
    bins: int
        how many the profiles were cut into, for labelling

    Returns
    -------
    `matplotlib.axes.Axes`
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6.4, 4.0))

    positions = np.arange(bins)
    width = 0.8 / max(len(profiles), 1)
    for index, (label, profile) in enumerate(profiles.items()):
        ax.bar(positions + index * width - 0.4 + width / 2, profile, width,
               label=label)

    ax.axhline(m_cut, color="0.4", linewidth=1, linestyle="--")
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{100 * i // bins}-{100 * (i + 1) // bins}%"
                        for i in range(bins)])
    ax.set_xlabel("fifth of the region's probability, least likely first")
    ax.set_ylabel("mean telescopes on a direction")
    # One series needs no legend, and a legend box over five bars is in the way.
    if len(profiles) > 1:
        ax.legend(frameon=False)
    return ax


def sky_bands(bands, source=None, centre_ra=None, ax=None, color="#2a78d6",
              title=None):
    """
    Credible bands of one localization, all-sky.

    Shaded by credible level rather than by raw probability density. Density
    spans orders of magnitude and plots as a thin bright thread with nothing
    around it; the bands answer the question actually being asked, which is
    where the source is likely to be and how confident that is.

    An all-sky Mollweide, because a region that runs a hundred degrees does not
    fit in anything smaller. For the close-up that keeps areas true, see
    `multiplicity_over_region`.

    Parameters
    ----------
    bands: dict of float to `divtel.region.SkyRegion`
        one region per credible level, in the sky's own frame -- not placed.
        Drawn outermost first so the darker inner bands land on top.
    source: `astropy.coordinates.SkyCoord`, optional
        a true position, marked if given
    centre_ra: `astropy.Quantity`, optional
        right ascension at the middle of the panel. Pass the same value across
        panels so several maps can be compared; the innermost band's mean
        otherwise.
    ax: `matplotlib.axes.Axes`, optional
        must have a Mollweide projection
    color: str
        the hue the bands are shaded in, light to dark
    title: str, optional

    Returns
    -------
    `matplotlib.axes.Axes`
    """
    from matplotlib.colors import to_rgb

    levels = sorted(bands)
    if ax is None:
        ax = plt.figure(figsize=(6.6, 4.4)).add_subplot(projection="mollweide")
    if centre_ra is None:
        centre_ra = bands[levels[0]].coord.ra.mean()

    base = to_rgb(color)
    # One hue, light to dark: the bands are a magnitude, not several identities.
    shades = dict(zip(levels,
                      [tuple(1 - (1 - channel) * weight for channel in base)
                       for weight in np.linspace(0.30, 1.0, len(levels))[::-1]],
                      strict=True))

    for depth, level in enumerate(reversed(levels)):
        coords = bands[level].coord
        ax.scatter(_wrap_ra(coords.ra, centre_ra), coords.dec.rad, s=1.4,
                   marker=".", linewidths=0, color=shades[level],
                   zorder=2 + depth, label=f"{level:.0%} credible")

    if source is not None:
        ax.plot(_wrap_ra(source.ra, centre_ra), source.dec.rad, marker="*",
                markersize=11, color="k", linestyle="none",
                markeredgecolor="w", markeredgewidth=0.7, zorder=10)

    ax.grid(True, linewidth=0.5)
    offsets = np.arange(-120, 121, 60)
    ax.set_xticks(np.radians(-offsets))
    ax.set_xticklabels([f"{(centre_ra.to_value(u.deg) + offset) % 360:.0f}$\\degree$"
                        for offset in offsets], fontsize=7)
    ax.tick_params(labelsize=7)
    if title:
        ax.set_title(title, pad=12)
    return ax


def _wrap_ra(ra, centre_ra):
    """Right ascension as a Mollweide x, with the panel centred on `centre_ra`."""
    # Mollweide runs right ascension right-to-left, as the sky is seen.
    return -(ra - centre_ra).wrap_at(180 * u.deg).rad

"""
Functions to define telescopes pointings
We use the same reference frame as simtel_array:
X is pointing North
Y is pointing West
Z is pointing upward
Az is taken from X (North) towards East, and between -180 and 180 degrees
Alt is taken from ground (towards Z) and between -90 and 90 degrees

East is -Y, so a pointing vector is (cos(alt)cos(az), -cos(alt)sin(az),
sin(alt)). This is astropy's azimuth convention too, which is what lets
`divtel.observation` hand alt/az straight to `Array.divergent_pointing` with no
conversion in between.
"""

import numpy as np
import astropy.units as u


def alt_az_to_vector(alt, az):
    """
    Compute a pointing vector coordinates (x,y,z) from an alt,az pointing direction

    Parameters
    ----------
    alt: float
        angle in rad
    az: float
        angle in rad

    Returns
    -------
    vector: `numpy.array`
        [x, y, z]
    """
    x = np.cos(alt.to(u.rad)) * np.cos(az.to(u.rad))
    y = -np.cos(alt.to(u.rad)) * np.sin(az.to(u.rad))
    z = np.sin(alt.to(u.rad))
    return np.array([x, y, z])


def as_altaz(vector):
    """
    A unit pointing vector back as an (alt, az) pair.

    The inverse of `alt_az_to_vector`.

    Parameters
    ----------
    vector: `numpy.ndarray`
        [x, y, z], unit length

    Returns
    -------
    (alt, az): tuple of `astropy.Quantity`, in degrees
    """
    alt = (np.arcsin(np.clip(vector[2], -1.0, 1.0)) * u.rad).to(u.deg)
    az = (np.arctan2(-vector[1], vector[0]) * u.rad).to(u.deg)
    return alt, az


def local_frame(alt, az):
    """
    Unit tangent vectors at an alt/az pointing, towards increasing az and alt.

    These are the directions an observer standing at that pointing and
    looking up would call "sideways" (towards increasing azimuth) and "up"
    (towards increasing altitude). They are exact, not approximations valid
    only near the pointing: alt/az is an orthogonal coordinate system on the
    sphere, so the two are always unit length and always perpendicular to
    each other and to `alt_az_to_vector(alt, az)`, at every pointing.

    Used to orient a small flat map -- a camera's field of view, or the sky
    map `Array.hyper_fov` builds -- around a pointing, so that moving along
    one returned vector reads as "more altitude" and the other as "more
    azimuth", whatever the pointing itself is.

    Parameters
    ----------
    alt, az: `astropy.Quantity`

    Returns
    -------
    (increasing_az, increasing_alt): tuple of `numpy.array`
        unit vectors [x, y, z]
    """
    alt = alt.to_value(u.rad)
    az = az.to_value(u.rad)
    increasing_az = np.array([-np.sin(az), -np.cos(az), 0.0])
    increasing_alt = np.array([-np.sin(alt) * np.cos(az),
                               np.sin(alt) * np.sin(az),
                               np.cos(alt)])
    return increasing_az, increasing_alt


def _norm_div(div, scale=100 * u.m):
    """
    Transformation function from div parameter to norm to compute the position of g_point

     Parameters
    ----------
    div: float
    scale: `astropy.Quantity`
        telescope distance from barycenter at which div = sin(divergence_angle)

    Returns
    -------
    `astropy.Quantity`
        distance, in metres
    """
    return scale/np.tan(np.arcsin(div))


def pointG_position(barycenter, div, alt_mean, az_mean):
    """
    Compute the position of g_point for the pointing

    Parameters
    ----------
    barycenter: `astropy.Quantity` or np.array([x,y,z])
        position of the barycenter of the array; a plain array is read as metres
    div: float
    alt_mean: `astropy.Quantity`
        mean pointing altitude in radians from which to diverge
    az_mean: `astropy.Quantity`
        mean pointing azimuth in radians from which to diverge

    Returns
    -------
    `astropy.Quantity` [Gx, Gy, Gz], in metres
    """
    barycenter = u.Quantity(barycenter, u.m)
    norm = _norm_div(div)
    g_x = barycenter[0] - norm * np.cos(alt_mean) * np.cos(az_mean)
    g_y = barycenter[1] + norm * np.cos(alt_mean) * np.sin(az_mean)
    g_z = barycenter[2] - norm * np.sin(alt_mean)
    return u.Quantity([g_x, g_y, g_z])


def tel_div_pointing(tel_position, g_point):
    """
    Divergent pointing to a point G.
    Update telescope pointing

    Parameters
    ----------
    tel_position: `astropy.Quantity` or np.array([x, y, z])
        telescope coordinates; a plain array is read as metres
    g_point: `astropy.Quantity` or numpy.array([Gx, Gy, Gz])

    Returns
    -------
    (alt, az): tuple of `astropy.Quantity`
        pointing direction, in radians
    """
    tel_position = u.Quantity(tel_position, u.m)
    g_point = u.Quantity(g_point, u.m)
    GT = np.sqrt(((tel_position - g_point) ** 2).sum())
    alt_tel = np.arcsin((tel_position[2] - g_point[2]) / GT)
    az_tel = np.arctan2(-(tel_position[1] - g_point[1]), (tel_position[0] - g_point[0]))
    return alt_tel.to(u.rad), az_tel.to(u.rad)


# The distance from the barycenter at which `div` is defined as the sine of the
# divergence angle -- the `scale` of `_norm_div`. Named here because it is part
# of the published meaning of div, not an implementation detail, and
# `div_for_half_angle` inverts it.
DIV_SCALE = 100 * u.m


def pointing_spread(array):
    """
    How far off the mean the array's telescopes are pointing.

    The observable side of the divergence parameter. `div` is a number between 0
    and 1 whose meaning is tied to a hundred metres of baseline; this is the
    angle on the sky it produces for the array actually in hand, which is what
    can be compared against the size of a target region.

    Parameters
    ----------
    array: `divtel.telescope.Array`
        already pointed

    Returns
    -------
    dict
        ``mean``, ``max`` and ``rms`` offset from the mean pointing direction,
        as `astropy.Quantity` angles
    """
    vectors = array.pointing_vectors
    mean = vectors.mean(axis=0)
    norm = np.linalg.norm(mean)
    mean = mean / norm if norm > 1e-9 else np.array([0.0, 0.0, 1.0])

    offsets = np.arccos(np.clip(vectors @ mean, -1.0, 1.0)) * u.rad

    return {
        "mean": offsets.mean().to(u.deg),
        "max": offsets.max().to(u.deg),
        "rms": np.sqrt((offsets**2).mean()).to(u.deg),
    }


def div_for_half_angle(array, half_angle, alt, az, tolerance=1e-4):
    """
    The divergence that spreads an array's pointings over a given half-angle.

    The direct answer to "the region is this big, how far do I have to diverge".
    A divergent pointing puts a point G a distance

        norm = 100 m / tan(arcsin(div))

    behind the array along the mean pointing, then aims every telescope along the
    line from G through itself. A telescope lying a perpendicular distance ``d``
    from the array's axis and a distance ``s`` along it therefore ends up

        arctan(d / (norm + s))

    off the mean pointing, and the array's spread is the largest of those. That
    is monotone decreasing in norm, so the norm putting the widest telescope
    exactly at `half_angle` is found by bisection, and

        div = sin(arctan(100 m / norm))

    puts G there.

    Both distances are measured against the pointing direction rather than
    against the ground, so an array seen edge-on at low elevation is correctly
    foreshortened and the same request costs more divergence than it would at the
    zenith. Keeping the along-axis term matters: dropping it, as a small-angle
    treatment does, understates the divergence needed by a third by the time the
    spread reaches 40 degrees, because the telescopes on the far side of the
    barycenter sit closer to G than the perpendicular distance alone suggests.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    half_angle: `astropy.Quantity`
        angle the outermost telescope should sit at, off the mean pointing
    alt, az: `astropy.Quantity`
        the mean pointing the array will diverge from
    tolerance: float
        convergence on div

    Returns
    -------
    float
        the div parameter, clipped into [0, 1]

    Notes
    -----
    Exact for the geometry, but it answers a question about the array's *widest*
    telescope, and it measures angles from the requested pointing while
    `pointing_spread` measures them from the mean of the pointings the array ends
    up with -- G lies behind the array, so the two differ slightly. Read the
    result back with `pointing_spread`. `solve_div` sidesteps all of this by
    searching on the covered weight itself.
    """
    axis = alt_az_to_vector(alt, az).astype(float)
    offsets = (array.positions_array - array.barycenter).to_value(u.m)

    along = offsets @ axis
    perpendicular = np.linalg.norm(offsets - np.outer(along, axis), axis=1)

    if perpendicular.max() <= 0:
        return 0.0

    target = half_angle.to_value(u.rad)

    def spread(norm):
        """Widest telescope offset, in radians, with G that far behind."""
        return np.arctan2(perpendicular, norm + along).max()

    # Bracket: G far away gives a parallel array, G at the barycenter gives a
    # spread of 90 degrees or more, so any reachable target lies between.
    low, high = float(np.abs(along).max()) + 1.0, 1e9
    if spread(low) < target:
        return 1.0

    while high - low > tolerance * low:
        middle = 0.5 * (low + high)
        if spread(middle) > target:
            low = middle
        else:
            high = middle

    return float(np.clip(np.sin(np.arctan(DIV_SCALE.to_value(u.m) / high)),
                         0.0, 1.0))


def div_for_multiplicity(array, alt, az, target=2.0, div_max=1.0, tolerance=1e-3):
    """
    The most an array can be spread before showers stop being seen twice.

    The ceiling on divergence, and it has nothing to do with the target. Spreading
    trades depth for width, and the trade has a hard floor: a shower seen by one
    telescope cannot be reconstructed stereoscopically, so once the mean
    multiplicity falls to two the array is spent. Past that point further
    divergence covers more sky with telescopes that can no longer do anything
    with what they see.

    This matters most for a small array. Thirteen telescopes asked to span sixty
    degrees end up with a mean multiplicity of one and cover nothing usefully,
    even though the geometry did exactly what was requested -- so a strategy that
    sizes divergence from the region alone will happily recommend a configuration
    that returns nothing.

    Mean multiplicity falls monotonically with divergence, so a bisection finds
    the crossing.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    alt, az: `astropy.Quantity`
        the pointing to evaluate at; the answer is mildly elevation-dependent,
        through the same foreshortening `div_for_half_angle` accounts for
    target: float
        multiplicity to hold, area-weighted over the sky the array sees. Two is
        the stereoscopic floor.
    div_max: float
        largest divergence considered
    tolerance: float
        convergence on div

    Returns
    -------
    float
        the divergence at which the mean multiplicity reaches `target`, or
        `div_max` if it never falls that far

    Notes
    -----
    Uses `divtel.telescope.Array.multiplicity_moments`, which averages over the
    sky the array covers rather than over any region -- the ceiling is a property
    of the configuration, not of what it happens to be pointed at.
    """
    def mean_multiplicity(div):
        array.divergent_pointing(div, alt, az)
        return array.multiplicity_moments()[0]

    if mean_multiplicity(div_max) >= target:
        return div_max

    low, high = 0.0, div_max
    while high - low > tolerance:
        middle = 0.5 * (low + high)
        if mean_multiplicity(middle) >= target:
            low = middle
        else:
            high = middle

    return low


def best_pointing(array, region, div=0.0, m_cut=2, rings=6, spokes=12, refine=6,
                  search_pixels=2000, seed=0, start=None):
    """
    Where to point, at a given divergence, to cover as much weight as possible.

    An array cannot always contain a region, and when it cannot, where it points
    decides how much it gets -- so the pointing is searched rather than assumed.
    Neither obvious guess survives contact with a real alert: for GW170817's
    first map the weighted centroid of the arc lies off the arc, and the axis of
    the smallest enclosing cone lies 37 degrees from the source.

    The candidates are a polar grid on the sky around the region's centroid,
    reaching out to the region's own extent plus one field of view. That bound is
    deliberate and it is not a performance compromise. Left unbounded, the search
    stops answering the question: a widely diverged array covers a ring rather
    than a disc, so aiming it sixty degrees off the region can land part of that
    ring squarely on it and score full marks, while every telescope not on the
    ring stares at empty sky. That configuration maximises the number this
    function returns and is worthless -- it is a handful of telescopes on the
    target and forty pointed nowhere, which is not divergent pointing but a badly
    chosen sub-array. Confining the aim to the region keeps "covered weight" a
    measure of how well the array is used, not of how cleverly it can be
    positioned to game a metric.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    region: `divtel.region.SkyRegion`
        placed over the array
    div: float
        divergence to search at; the answer depends on it, since a widely spread
        configuration wants aiming differently from a tight one
    m_cut: int
    rings, spokes: int
        the polar grid: `rings` distances out to the array's reach, `spokes`
        directions around each
    refine: int
        rounds of local search after the grid, each halving the step
    search_pixels: int
        directions the search itself scores against. The grid is swept hundreds
        of times over, and ranking candidates does not need every one -- the
        winner is rescored against the whole region before being returned, so the
        number that comes back is exact whatever this is.
    seed: int
        for the subsample, so a rerun gives the same answer
    start: (`astropy.Quantity`, `astropy.Quantity`), optional
        an extra (alt, az) candidate to seed the search with. The best pointing
        moves smoothly with divergence, so handing a scan's previous answer back
        in costs one evaluation and stops neighbouring points on the curve landing
        in different local optima -- which shows up as a curve that jumps about
        rather than as an obviously wrong answer.

    Returns
    -------
    dict
        ``alt``, ``az``, ``covered``; the array is left pointed at the best one
    """
    directions, weights = region.directions, region.weights

    centre = np.average(directions, axis=0, weights=weights)
    centre /= np.linalg.norm(centre)

    if len(directions) > search_pixels:
        chosen = np.random.default_rng(seed).choice(
            len(directions), search_pixels, replace=False)
        sample, sample_weights = directions[chosen], weights[chosen]
    else:
        sample, sample_weights = directions, weights

    cos_radii = np.cos([t.fov_radius.to_value(u.rad) for t in array.telescopes])

    def aim(vector):
        """Point the array along a vector, and say how much it then covers."""
        vector = vector / np.linalg.norm(vector)
        alt, az = as_altaz(vector)
        array.divergent_pointing(div, alt, az)
        multiplicity = ((sample @ array.pointing_vectors.T)
                        >= cos_radii).sum(axis=1)
        covered = float(sample_weights[multiplicity >= m_cut].sum()
                        / sample_weights.sum())
        return covered, vector, alt, az

    # The aim may wander over the region and one field of view beyond it, and no
    # further -- see the note above on why this is bounded.
    reach = (np.arccos(np.clip((directions @ centre).min(), -1.0, 1.0))
             + max(t.fov_radius for t in array.telescopes).to_value(u.rad))

    east, north = local_frame(*as_altaz(centre))

    best = aim(centre)
    if start is not None:
        trial = aim(alt_az_to_vector(*start).astype(float))
        if trial[0] > best[0]:
            best = trial

    for radius in np.linspace(reach / rings, reach, rings):
        for angle in np.linspace(0, 2 * np.pi, spokes, endpoint=False):
            trial = aim(centre
                        + radius * (np.cos(angle) * east + np.sin(angle) * north))
            if trial[0] > best[0]:
                best = trial

    step = reach / rings
    for _ in range(refine):
        improved = False
        for direction in (east, -east, north, -north,
                          east + north, east - north, -east + north, -east - north):
            trial = aim(best[1] + step * direction / np.linalg.norm(direction))
            if trial[0] > best[0]:
                best, improved = trial, True
        if not improved:
            step *= 0.5

    _, _, alt, az = best
    array.divergent_pointing(div, alt, az)
    return {"alt": alt, "az": az, "covered": region.covered(array, m_cut)}


def div_scan(array, region, divs, alt=None, az=None, m_cuts=(1, 2),
             hyper_fov=True, optimize_pointing=False):
    """
    Coverage, field of view and multiplicity across a range of divergences.

    The curve the divergence trade turns on: as div grows the array covers more
    of the region, and fewer telescopes see any part of it.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    region: `divtel.region.SkyRegion`
        placed over the array
    divs: iterable of float
        divergences to try
    alt, az: `astropy.Quantity`, optional
        mean pointing to diverge from. Defaults to the region's centroid, which
        is where an observer with no better idea would point.
    m_cuts: iterable of int
        multiplicity cuts to report coverage for
    hyper_fov: bool
        also compute the array's total covered solid angle and multiplicity
        moments. Accurate but slow -- it is the polygon arrangement in
        `divtel.telescope.Array.hyper_fov` -- so it can be turned off for a fine
        scan.
    optimize_pointing: bool
        re-choose the pointing at each divergence with `best_pointing`, rather
        than holding `alt` and `az` fixed. What a real observer would do, and it
        matters for a region the array cannot contain: the best place to point a
        tightly-packed array is not the best place to point a spread one.

    Returns
    -------
    list of dict
        one per divergence, with ``div``, ``alt``, ``az``, ``spread``,
        ``covered`` (a dict keyed by m_cut), and when asked ``hyper_fov_area``,
        ``stereo_area``, ``multiplicity_mean`` and ``multiplicity_std``
    """
    m_cuts = tuple(m_cuts)
    if alt is None or az is None:
        alt, az = as_altaz(region.centroid)

    weights = region.weights
    rows = []
    previous = None

    for div in divs:
        if optimize_pointing:
            chosen = best_pointing(array, region, float(div), max(m_cuts),
                                   start=previous)
            pointing = (chosen["alt"], chosen["az"])
            previous = pointing
        else:
            array.divergent_pointing(float(div), alt, az)
            pointing = (alt, az)

        multiplicity = region.multiplicity(array)

        row = {
            "div": float(div),
            "alt": pointing[0],
            "az": pointing[1],
            "spread": pointing_spread(array),
            "covered": {m: float(weights[multiplicity >= m].sum() / weights.sum())
                        for m in m_cuts},
        }

        if hyper_fov:
            area, patches = array.hyper_fov(min_telescopes=1)
            stereo, _ = array.hyper_fov(min_telescopes=2)
            mean, variance = array.multiplicity_moments(patches)
            row["hyper_fov_area"] = area
            row["stereo_area"] = stereo
            row["multiplicity_mean"] = mean
            row["multiplicity_std"] = float(np.sqrt(variance))

        rows.append(row)

    return rows


def spread_scan(array, divs, alt, az):
    """
    How divergence trades depth for width, with no region involved.

    `div_scan` folds a target in at every step, because covered weight is what a
    strategy is chasing. But two of the numbers it reports along the way -- mean
    multiplicity and the hyper field of view's area -- do not need a target to
    mean something. They describe the array and the divergence alone.

    The two move in lockstep, and by a fixed amount. Multiplicity is a sum of
    camera discs, so ``<m> = total camera solid angle / hyper_fov area``
    regardless of how the array is spread -- that identity is what
    `Array.multiplicity_moments` computes. Dividing the hyper field of view by
    the array's mean single-camera area instead of by its raw solid angle turns
    that into ``<m> * n_fov = n_telescopes``, a constant: as the array spreads,
    ``n_fov`` grows from one camera's worth of sky towards `n_telescopes`
    cameras' worth exactly as fast as `<m>` falls from `n_telescopes` towards
    one. Neither curve carries information the other does not; the point of
    plotting both is that "one camera wide, seen fifty times over" and "fifty
    cameras wide, seen once" are the two ends of the same trade, in units an
    observer can picture.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    divs: iterable of float
        divergences to try
    alt, az: `astropy.Quantity`
        mean pointing to diverge from; the hyper field of view is mildly
        elevation-dependent through the same foreshortening `div_for_half_angle`
        accounts for, so pick one pointing and hold it for the whole scan

    Returns
    -------
    list of dict
        one per divergence, with ``div``, ``spread``, ``multiplicity_mean``,
        ``multiplicity_std``, ``hyper_fov_area`` (`astropy.Quantity`, at
        min_telescopes 1) and ``n_fov`` -- the hyper field of view measured in
        cameras, i.e. ``hyper_fov_area`` divided by the array's mean
        single-camera solid angle
    """
    radii = np.array([t.fov_radius.to_value(u.rad) for t in array.telescopes])
    mean_camera_area = float(np.mean(2 * np.pi * (1 - np.cos(radii)))) * u.sr

    rows = []
    for div in divs:
        array.divergent_pointing(float(div), alt, az)
        area, patches = array.hyper_fov(min_telescopes=1)
        mean, variance = array.multiplicity_moments(patches)
        rows.append({
            "div": float(div),
            "spread": pointing_spread(array),
            "multiplicity_mean": mean,
            "multiplicity_std": float(np.sqrt(variance)),
            "hyper_fov_area": area,
            "n_fov": float((area / mean_camera_area)
                           .to_value(u.dimensionless_unscaled)),
        })

    return rows


def solve_div(array, region, alt=None, az=None, target=0.9, m_cut=2, coarse=41,
              tolerance=1e-3, optimize_pointing=False):
    """
    The least divergence that covers a given fraction of the region.

    Searched rather than derived. A coarse sweep of the whole range first, then
    bisection inside the bracket it finds: coverage rises with divergence at
    first and can fall again once the fields of view pull apart and the
    multiplicity cut starts biting, so bisecting from the outset could converge
    onto the wrong side of a maximum.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    region: `divtel.region.SkyRegion`
        placed over the array
    alt, az: `astropy.Quantity`, optional
        mean pointing to diverge from; the region's centroid by default
    target: float
        fraction of the region's weight to cover
    m_cut: int
    coarse: int
        points in the initial sweep of div over [0, 1]
    tolerance: float
        stop when the bracket on div is this narrow
    optimize_pointing: bool
        re-aim at every divergence with `best_pointing`

    Returns
    -------
    dict
        ``div`` and ``covered`` at the solution, ``spread`` there, and
        ``best_div``/``best_covered`` -- the most the array can cover at all.
        ``div`` is None when the target is out of reach, and then the best
        entries say how close it gets and where.
    """
    if alt is None or az is None:
        alt, az = as_altaz(region.centroid)

    def covered(div):
        if optimize_pointing:
            return best_pointing(array, region, float(div), m_cut)["covered"]
        array.divergent_pointing(float(div), alt, az)
        return region.covered(array, m_cut)

    grid = np.linspace(0.0, 1.0, coarse)
    values = np.array([covered(div) for div in grid])

    best = int(np.argmax(values))
    result = {"best_div": float(grid[best]), "best_covered": float(values[best])}

    reached = np.flatnonzero(values >= target)
    if len(reached) == 0:
        covered(grid[best])
        result.update({"div": None, "covered": float(values[best]),
                       "spread": pointing_spread(array),
                       "alt": array.mean_pointing[0], "az": array.mean_pointing[1]})
        return result

    first = reached[0]
    if first == 0:
        high = 0.0
    else:
        low, high = grid[first - 1], grid[first]
        while high - low > tolerance:
            middle = 0.5 * (low + high)
            if covered(middle) >= target:
                high = middle
            else:
                low = middle

    result.update({"div": float(high), "covered": covered(high),
                   "spread": pointing_spread(array),
                   "alt": array.mean_pointing[0], "az": array.mean_pointing[1]})
    return result

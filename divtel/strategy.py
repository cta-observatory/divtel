"""
Covering a weighted region of sky: tiling, sub-arrays, and shaped pointing.

Divergent pointing spreads an array with one number, and one number buys a
symmetric configuration -- a disc, or a ring of discs. Real targets are neither.
A two-detector gravitational-wave localization is an arc a hundred degrees long
and a few degrees wide with the probability piled unevenly along it, and no
value of ``div`` makes fifty telescopes look like that. This module holds the
alternatives, all of them working on a `divtel.region.SkyRegion`.

`tile_region`
    Keep the array together and visit the region in pieces. Nothing is wasted
    and every shower is seen by every telescope; what is spent is time, since
    the window has to be divided between the pointings.
`point_subarrays`, `point_subarrays_by_type`
    Split the array into groups and point each group conventionally at a
    different part of the region, all at once. Every telescope is on target,
    nothing waits its turn, and each group is an ordinary small array.
`shaped_pointing`, `shaped_pointing_by_type`
    Throw the groups away too and point every telescope on its own, so that the
    array's depth follows the weight rather than merely covering it.
`weighted_split`
    The middle road: keep the sub-array framework but deal the *sizes* in
    proportion to the weight each group catches. Most of shaped pointing's
    gradient, without a solver and without a table of fifty directions.

The rule behind shaped pointing
-------------------------------
Coverage is not the objective. What the array delivers at a direction is its
*multiplicity* there -- how many cameras have it in frame -- and that is a
continuous resource which can be spent unevenly. Two demands:

    1. see all of the region, not part of it
    2. look hardest where the weight is

Both fall out of one choice. Maximise

    J = sum_k  w_k  log( max(m_k - m_cut + 1, 0) + eps )

over the pointings, where `w` is the weight of direction `k` and `m` its
multiplicity. The logarithm is what makes this work, and it is not a fudge:

* Its optimum under a fixed budget is **m proportional to w**. Maximising
  `sum w log m` subject to `sum m = B` gives `m_k = B w_k` and nothing else.
  Demand 2 is not imposed, it is what a logarithm returns.
* It diverges downward at zero, so leaving any weight unseen costs infinitely
  much and demand 1 comes free. `eps` is the finite stand-in, and turning it up
  is how you buy depth on the core at the price of the tails.
* Subtracting `m_cut - 1` inside means a direction seen by one telescope is
  worth exactly what an unseen one is worth. A shower nobody triangulates is not
  a detection, and an optimiser told otherwise will spread the array into a thin
  single-telescope sheet that scores beautifully and reconstructs nothing. This
  is the term that stops it -- the same degenerate optimum divergence falls into
  past its multiplicity ceiling, reached from the other side.

Where the budget runs out
-------------------------
Demand 1 is not always satisfiable. The array has a fixed amount of camera to
spend -- `reach` measures it, as the multiplicity the array could hold over the
region if none of its field of view spilled off the edge -- and when that is
below `m_cut` no arrangement covers the whole region stereoscopically. Then the
right answer is to give part of it up, and to give up the *least* weighted part.
`target_field` does that explicitly, by water-filling.

That zero is a decision an observer has to make anyway. Writing it down beats
discovering it afterwards in a configuration that covered everything once.

What it costs
-------------
The configuration is no longer one number. It is a table of alt/az, one row per
telescope, which is what `divtel.telescope.Array.export_cfg` writes. And the
pointings are optimised against one weight map at one moment, so a map that is
superseded -- as gravitational-wave maps always are -- invalidates them in a way
a divergence value does not.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np

from .region import as_altaz

__all__ = [
    "stereo_radius",
    "tile_region",
    "split_into_subarrays",
    "point_subarrays",
    "point_subarrays_by_type",
    "weighted_split",
    "reach",
    "target_field",
    "shaped_pointing",
    "shaped_pointing_by_type",
    "describe",
    "camera_blur",
    "multiplicity_by_probability",
]


def stereo_radius(array, m_cut=2):
    """
    The angular radius a normally-pointed array covers deeply enough to use.

    Pointed normally every telescope looks the same way, so a direction is seen
    by whichever cameras are wide enough to reach it. Sorting the cameras by
    radius, a direction is seen by `m_cut` or more of them exactly when it lies
    within the `m_cut`-th widest.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    m_cut: int

    Returns
    -------
    `astropy.Quantity`
        angular radius

    Raises
    ------
    ValueError
        if the array has fewer than `m_cut` telescopes
    """
    radii = sorted((t.fov_radius for t in array.telescopes), reverse=True)
    if len(radii) < m_cut:
        raise ValueError(
            f"an array of {len(radii)} telescopes cannot see anything "
            f"{m_cut} times over"
        )
    return radii[m_cut - 1]


def reach(array, region_area):
    """
    The multiplicity the array could hold over a region, if nothing spilled.

    Total camera solid angle divided by the region's, so it is the mean
    multiplicity of a perfect arrangement: every camera entirely on the region
    and the overlaps distributed exactly as asked. Nothing reaches it -- a disc
    laid on the edge of a region hangs over the side, and the loss grows with
    the camera -- but it is an upper bound that costs one division to compute,
    and it is the number to look at before anything else.

    Below 1 the region cannot be covered at all. Below `m_cut` it cannot be
    covered stereoscopically, and part of it has to be abandoned whatever the
    strategy. `shaped_pointing` reports the fraction actually achieved as
    ``efficiency``.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    region_area: `astropy.Quantity`
        solid angle of the region, as `divtel.region.SkyRegion.area` gives it

    Returns
    -------
    float
        mean multiplicity available, dimensionless
    """
    radii = np.array([t.fov_radius.to_value(u.rad) for t in array.telescopes])
    camera = (2 * np.pi * (1 - np.cos(radii))).sum()
    return float(camera / region_area.to_value(u.sr))


def tile_region(array, region, pointings=10, m_cut=2, candidates=4000, seed=0):
    """
    Cover a region with a sequence of ordinary pointings.

    The alternative to spreading the array: leave it pointed normally and move
    it, taking the region in pieces. Each pointing puts every telescope on the
    same patch, so nothing is wasted on empty sky and the multiplicity stays at
    its maximum. What is spent instead is time, since the observing window has
    to be divided between the pointings.

    Tiles are chosen greedily. The first goes where it catches the most weight,
    and each one after it where it catches the most of what is still uncovered.
    Greedy disc covering is not optimal, but it is within a small constant factor
    of optimal and it matches how an observer would actually work down a region.

    Parameters
    ----------
    array: `divtel.telescope.Array`
        used only for its camera sizes; it is not re-pointed
    region: `divtel.region.SkyRegion`
        placed over the array
    pointings: int
        how many tiles to plan
    m_cut: int
        telescopes that must see a direction for it to count
    candidates: int
        tile centres considered, drawn from the region's own directions
    seed: int
        for the subsample, so a rerun gives the same tiling

    Returns
    -------
    list of dict
        one per tile in order, with ``alt``, ``az``, ``new`` (weight this tile
        adds, as a fraction of the region) and ``cumulative``
    """
    directions, weights = region.directions, region.weights
    total = weights.sum()
    cos_radius = np.cos(stereo_radius(array, m_cut).to_value(u.rad))

    index = np.arange(len(directions))
    if len(index) > candidates:
        index = np.random.default_rng(seed).choice(len(index), candidates,
                                                   replace=False)

    # Which directions each candidate tile would cover. One matrix, reused by
    # every round of the greedy loop, so the cost is paid once.
    covers = (directions[index] @ directions.T) >= cos_radius

    covered = np.zeros(len(directions), dtype=bool)
    tiles = []
    for _ in range(pointings):
        gain = (covers & ~covered) @ weights
        best = int(np.argmax(gain))
        if gain[best] <= 0:
            break

        covered |= covers[best]
        alt, az = as_altaz(directions[index[best]])
        tiles.append({
            "alt": alt, "az": az,
            "new": float(gain[best] / total),
            "cumulative": float(weights[covered].sum() / total),
        })

    return tiles


def split_into_subarrays(array, count):
    """
    Deal an array into sub-arrays of even size and even reach.

    Called on one telescope type at a time by `point_subarrays_by_type`, which
    is the normal path: an MST and an SST are different instruments with
    different fields of view, and mixing them in one sub-array was the earlier
    version of this idea. A sub-array that happened to collect only the narrower
    type reached less far than one that collected the wider, and the sub-array
    with the narrowest reach set what the whole plan could do.

    Nothing here actually requires a single type, though. Telescopes are sorted
    by field-of-view radius and dealt round-robin, the way cards are dealt, so
    every sub-array ends up the same size and the same mix of whatever radii went
    in. Feed it a mixed array and that mix balances the reach; feed it one type
    and the sorting has nothing left to do.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    count: int
        how many sub-arrays to make

    Returns
    -------
    dict of str to `divtel.telescope.Array`
        named ``"1"``, ``"2"`` and so on. The sub-arrays hold the same
        `Telescope` objects as the parent, so pointing one points the parent's
        telescopes too, and multiplicity can be measured on the whole array.

    Raises
    ------
    ValueError
        if there are not enough telescopes to fill every sub-array
    """
    if count < 1:
        raise ValueError(f"need at least one sub-array, got {count}")
    if count > len(array.telescopes):
        raise ValueError(
            f"cannot deal {len(array.telescopes)} telescopes into {count} "
            "sub-arrays without leaving some empty"
        )

    ordered = sorted(array.telescopes, key=lambda t: -t.fov_radius)
    groups = {str(n + 1): [] for n in range(count)}
    for position, telescope in enumerate(ordered):
        groups[str(position % count + 1)].append(telescope.id)

    return array.group_by(groups)


def _place_groups(array, region, groups, m_cut, candidates, seed):
    """
    Point each sub-array greedily at the weight the others have not taken.

    Widest reach first: the sub-array that can cover the most should choose from
    the whole region rather than from what narrower ones left behind.
    """
    directions, weights = region.directions, region.weights
    total = weights.sum()

    order = sorted(groups, key=lambda name: -stereo_radius(groups[name], m_cut)
                   if len(groups[name].telescopes) >= m_cut else np.inf)

    index = np.arange(len(directions))
    if len(index) > candidates:
        index = np.random.default_rng(seed).choice(len(index), candidates,
                                                   replace=False)
    overlap = directions[index] @ directions.T

    covered = np.zeros(len(directions), dtype=bool)
    pointings = []
    for name in order:
        group = groups[name]
        if len(group.telescopes) < m_cut:
            # Too few telescopes to see anything stereoscopically. Park it on
            # the densest uncovered spot rather than pretending it contributes.
            covers = np.zeros_like(overlap, dtype=bool)
        else:
            covers = overlap >= np.cos(
                stereo_radius(group, m_cut).to_value(u.rad))

        gain = (covers & ~covered) @ weights
        best = int(np.argmax(gain))
        covered |= covers[best]

        alt, az = as_altaz(directions[index[best]])
        group.divergent_pointing(0.0, alt, az)
        pointings.append({"subarray": name, "telescopes": len(group.telescopes),
                          "alt": alt, "az": az,
                          "new": float(gain[best] / total)})

    return pointings


def _summarise_groups(array, region, groups, pointings, m_cut):
    """What a set of pointed sub-arrays delivers, measured on the whole array."""
    directions, weights = region.directions, region.weights

    # Measured on the whole array, so a direction two sub-arrays both happen to
    # see is counted once and with the multiplicity it really gets.
    multiplicity = region.multiplicity(array)
    seen = multiplicity >= m_cut

    cos_radii = np.cos([t.fov_radius.to_value(u.rad) for t in array.telescopes])
    on_target = int(((directions @ array.pointing_vectors.T) >= cos_radii)
                    .any(axis=0).sum())

    return {
        "count": len(groups),
        "covered": float(weights[seen].sum() / weights.sum()),
        "multiplicity": (float(np.average(multiplicity[seen], weights=weights[seen]))
                         if seen.any() else 0.0),
        "on_target": on_target,
        "sizes": [len(group.telescopes) for group in groups.values()],
        "pointings": pointings,
    }


def point_subarrays(array, region, count, m_cut=2, candidates=4000, seed=0):
    """
    Split the array and aim each piece at a different part of the region.

    The configuration an observer would actually reach for. Divergent pointing
    spreads one array thin over a patch that is mostly empty sky; tiling covers
    the region properly but only one piece at a time. Sub-arrays do both at once:
    each piece is a small ordinary array staring at one part of the region, all
    of them at the same time, for the whole window.

    What it costs is telescopes per shower. Three sub-arrays of seventeen see
    every shower seventeen times instead of fifty-one, which is still far above
    the two a stereoscopic reconstruction needs.

    Sub-arrays are placed greedily, the one with the widest reach first, each
    going where it catches the most weight not already covered. Reach is computed
    per sub-array rather than assumed: dealt thinly enough, a sub-array runs out
    of wide cameras and its reach collapses to a narrow one, which is a real
    limit on how far this can be pushed.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    region: `divtel.region.SkyRegion`
        placed over the array
    count: int
        how many sub-arrays to make
    m_cut: int
        telescopes within one sub-array that must see a direction
    candidates: int
        pointings considered, drawn from the region's own directions
    seed: int

    Returns
    -------
    dict
        ``covered`` (fraction of the region seen by `m_cut` telescopes or more,
        counting the whole array), ``multiplicity`` (weighted mean over the
        covered part), ``on_target`` (telescopes with the region in frame),
        ``sizes`` and ``pointings``. The array is left pointed.
    """
    groups = split_into_subarrays(array, count)
    pointings = _place_groups(array, region, groups, m_cut, candidates, seed)
    ordered = {name: groups[name] for name in sorted(groups, key=int)}
    return _summarise_groups(array, region, ordered, pointings, m_cut)


def point_subarrays_by_type(array, region, target=0.9, m_cut=2, candidates=4000,
                            seed=0):
    """
    Split each telescope type on its own, and let each cover the region by itself.

    `point_subarrays` deals telescopes into balanced sub-arrays precisely to stop
    a mixed sub-array reaching only as far as its narrowest telescope. The
    alternative is not to balance the mix but to drop it: two telescope types are
    different instruments, and there is no reason the number of pieces that suits
    one has to suit the other. Fourteen MSTs might split three ways to reach the
    region; thirty-seven SSTs, ten times the camera, might need only two. Solving
    that as one array with one sub-array count forces a compromise neither type
    asked for.

    So each type is handed the whole problem separately: split by field-of-view
    radius, and within each type scan sub-array counts from one upward until
    `target` coverage is reached or the type runs out of telescopes to add. A
    type with fewer than `m_cut` telescopes cannot see anything stereoscopically
    on its own and is left pointed nowhere in particular.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    region: `divtel.region.SkyRegion`
    target: float
        coverage each type's scan aims for before it stops adding sub-arrays
    m_cut: int
    candidates: int
    seed: int

    Returns
    -------
    dict
        ``types``, one entry per field-of-view radius holding that type's chosen
        `point_subarrays` result (or ``None`` if the type has fewer than `m_cut`
        telescopes); ``covered`` and ``multiplicity``, measured on the whole
        array once every type is pointed, so a direction two types both reach
        counts once at its true combined multiplicity rather than twice.
    """
    weights = region.weights

    types = {}
    for name, group in array.group_by("fov_radius").items():
        if len(group.telescopes) < m_cut:
            types[name] = None
            continue

        rows = [point_subarrays(group, region, count, m_cut, candidates, seed)
                for count in range(1, len(group.telescopes) // m_cut + 1)]
        chosen = next((r for r in rows if r["covered"] >= target), rows[-1])
        # The scan leaves the group pointed at its last count, not the chosen
        # one -- repoint it before moving to the next type.
        types[name] = point_subarrays(group, region, chosen["count"], m_cut,
                                      candidates, seed)

    multiplicity = region.multiplicity(array)
    seen = multiplicity >= m_cut

    return {
        "types": types,
        "covered": float(weights[seen].sum() / weights.sum()),
        "multiplicity": (float(np.average(multiplicity[seen], weights=weights[seen]))
                         if seen.any() else 0.0),
    }


def weighted_split(array, region, count, m_cut=2, candidates=4000, seed=0):
    """
    Sub-arrays whose *sizes* follow the weight each of them catches.

    `split_into_subarrays` deals round-robin, so every group holds the same
    number of telescopes; the multiplicity inside a group's footprint is that
    number, and the depth profile is therefore flat by construction rather than
    by any property of tiling. A gradient needs unequal groups.

    This builds one. An even split is placed first, to fix where the groups go
    and how much weight each catches; the telescopes are then re-dealt so that
    each group's size is proportional to its catch, subject to every group
    keeping at least `m_cut`; and the groups are re-pointed at the same places.

    Holding the placement fixed means the answer is a lower bound on what unequal
    splits can do -- a joint optimisation over placement and size would do at
    least as well. Even so it recovers most of `shaped_pointing`'s gradient
    without a solver and without leaving the sub-array framework, which is a
    change to a scheduling heuristic rather than a new method.

    Parameters
    ----------
    array: `divtel.telescope.Array`
    region: `divtel.region.SkyRegion`
    count: int
        how many sub-arrays to make
    m_cut: int
    candidates: int
    seed: int

    Returns
    -------
    dict
        as `point_subarrays`, with ``sizes`` now unequal

    Raises
    ------
    ValueError
        if the array cannot give every group its `m_cut` floor
    """
    if count * m_cut > len(array.telescopes):
        raise ValueError(
            f"{count} sub-arrays each holding at least {m_cut} telescopes needs "
            f"{count * m_cut}, and the array has {len(array.telescopes)}"
        )

    even = split_into_subarrays(array, count)
    placement = _place_groups(array, region, even, m_cut, candidates, seed)
    caught = np.array([row["new"] for row in placement], dtype=float)
    names = [row["subarray"] for row in placement]

    # Largest-remainder apportionment above the floor, so the sizes are integers
    # that sum to the array exactly rather than to it plus rounding.
    surplus = len(array.telescopes) - count * m_cut
    share = caught / caught.sum() if caught.sum() > 0 else np.full(count, 1 / count)
    exact = share * surplus
    sizes = np.floor(exact).astype(int)
    for position in np.argsort(-(exact - sizes))[:surplus - int(sizes.sum())]:
        sizes[position] += 1
    sizes += m_cut

    # Re-deal widest cameras to the hungriest groups, so a group asked to be deep
    # is not handed only the narrow telescopes the previous group left over.
    ordered = sorted(array.telescopes, key=lambda t: -t.fov_radius)
    groups, cursor = {}, 0
    for name, size in zip(names, sizes, strict=True):
        groups[name] = [t.id for t in ordered[cursor:cursor + int(size)]]
        cursor += int(size)
    groups = array.group_by(groups)

    for row in placement:
        groups[row["subarray"]].divergent_pointing(0.0, row["alt"], row["az"])
        row["telescopes"] = len(groups[row["subarray"]].telescopes)

    ordered_groups = {name: groups[name] for name in sorted(groups, key=int)}
    return _summarise_groups(array, region, ordered_groups, placement, m_cut)


def target_field(weights, budget, beta=1.0, m_cut=2):
    """
    The multiplicity the strategy asks for at each direction.

    Water-filling. The most weighted directions are served first, each at
    `m_cut`, until the budget is spent; whatever remains is poured on top in
    proportion to weight. Directions never reached are written as zero, and that
    zero is the strategy saying which part of the region it is giving up on.

    With a generous budget the floor is a rounding detail and the field is
    essentially `budget * weights`: the proportional allocation the logarithmic
    objective wants. With a tight one the floor is everything and the field is
    flat over whatever fraction can be seen in stereo at all.

    Parameters
    ----------
    weights: `numpy.ndarray`
        weight per direction; need not be normalised
    budget: float
        total multiplicity to spend, in per-direction units -- `reach` times the
        number of directions. `shaped_pointing` re-estimates this from what it
        achieves rather than trusting the geometric value.
    beta: float
        how sharply to follow the weight once the floor is paid. 1 tracks it,
        which is the point of this module. 0 spreads the surplus evenly and asks
        only for coverage. Above 1 concentrates harder on the core than the
        weight warrants, which is the right call only if sensitivity rises faster
        than linearly with multiplicity.
    m_cut: int
        multiplicity below which a direction is worth nothing

    Returns
    -------
    `numpy.ndarray`
        target multiplicity per direction, summing to `budget` unless the floor
        alone exhausts it
    """
    weights = np.asarray(weights, dtype=float)
    weights = weights / weights.sum()
    n_pixels = len(weights)

    floor = max(float(m_cut), 0.0)
    if floor <= 0:
        shape = weights ** beta if beta > 0 else np.ones(n_pixels)
        return shape / shape.sum() * budget

    # How many directions the budget can hold at the floor, taken most weighted
    # first. Fewer than all of them is the interesting case and the honest one.
    served = int(min(n_pixels, max(1, budget // floor)))
    order = np.argsort(-weights)[:served]

    field = np.zeros(n_pixels)
    field[order] = floor
    surplus = budget - served * floor
    if surplus > 0:
        shape = weights[order] ** beta if beta > 0 else np.ones(served)
        field[order] += surplus * shape / shape.sum()
    return field


def _utility(multiplicity, m_cut, eps):
    """Value of a direction seen this many times. Flat and worthless below `m_cut`."""
    return np.log(np.maximum(multiplicity - (m_cut - 1), 0.0) + eps)


def _point(array, vectors):
    """Aim every telescope along its own vector."""
    for telescope, vector in zip(array.telescopes, vectors, strict=True):
        telescope.point_to_altaz(*as_altaz(vector))


def shaped_pointing(array, region, region_area=None, beta=1.0, m_cut=2, eps=0.1,
                    candidates=1500, search_pixels=2500, sweeps=12, rounds=3,
                    warm_starts=(2, 3, 4, 6, 8, 10, 14), seed=0):
    """
    Point each telescope where it does the most good, and leave it there.

    The solver behind the rule in this module's docstring. It maximises
    ``sum w log(max(m - m_cut + 1, 0) + eps)`` over the pointings by coordinate
    ascent: take one telescope out of the arrangement, work out what the sky
    looks like without it, and put it back wherever it now helps most. Sweep
    until nothing wants to move.

    That inner step is one matrix-vector product. Precompute which directions
    each candidate pointing would cover -- there are only two camera sizes at
    either CTAO site, so two matrices -- and the marginal gain of every candidate
    for one telescope is that matrix against the per-direction value of one more
    telescope. Fifty telescopes times a dozen sweeps is a few thousand of them.

    Every move is scored on the exact objective and only taken if it improves, so
    the sweeps climb monotonically and stop. Where they stop depends on where
    they started, which is why several starts are tried: a greedy build, and the
    sub-array configurations from `point_subarrays`, which are good arrangements
    arrived at by a different argument and hard for a local search to find on its
    own.

    The budget is not known in advance. `reach` assumes no camera hangs off the
    edge of the region and every real arrangement loses something that way, by a
    factor of two on a compact region and rather more on a thin arc. Trusting the
    geometric value asks for a floor the array cannot hold, and the answer is a
    region covered once over rather than half of it covered properly. So the
    budget is calibrated: solve, measure the multiplicity that actually landed,
    rebuild the target from that, solve again.

    Parameters
    ----------
    array: `divtel.telescope.Array`
        left pointed at the answer
    region: `divtel.region.SkyRegion`
        placed over the array
    region_area: `astropy.Quantity`, optional
        solid angle of that region, for the initial budget. Taken from the
        region itself when it knows its own pixel area.
    beta: float
        how hard to follow the weight; see `target_field`
    m_cut: int
        telescopes that must see a direction for it to be worth anything
    eps: float
        value of an unseen direction, as ``log(eps)``. Small insists on covering
        the tails before deepening the core; large trades the tails away. The
        default sits nearer the coverage end.
    candidates: int
        pointings considered per telescope, drawn from the region's own
        directions. The candidate set is the resolution limit on where a
        telescope can be put, and at HEALPix nside 256 the pixels are 0.2 degrees
        apart against cameras of several degrees, so this is not the binding
        approximation.
    search_pixels: int
        directions the sweeps score against. The result is rescored against the
        whole region before it is returned, so the numbers that come back are
        exact whatever this is.
    sweeps: int
        most passes over the telescopes per start
    rounds: int
        budget calibrations
    warm_starts: tuple of int
        sub-array counts to try as starting arrangements, alongside the greedy
        build. Cheap and worth having: a split array is a good arrangement
        arrived at by a different argument, and in the budget-limited case --
        where there is not enough camera to cover the region twice over and the
        answer is a patchwork rather than a gradient -- it is the start the
        sweeps most often end up improving on. Empty to skip them, which is
        faster and sometimes worse.
    seed: int
        for the subsamples, so a rerun gives the same answer

    Returns
    -------
    dict
        ``budget`` (mean multiplicity the arrangement puts on the region, the
        calibrated figure), ``reach`` (the geometric bound), ``efficiency``
        (their ratio), ``pointings`` (list of alt/az, one per telescope, in array
        order) and ``objective``. The metrics live in `describe`.
    """
    if region_area is None:
        region_area = region.area

    directions = region.directions
    weights = np.asarray(region.weights, dtype=float)
    weights = weights / weights.sum()
    n_pixels = len(directions)
    rng = np.random.default_rng(seed)

    radii = np.array([t.fov_radius.to_value(u.rad) for t in array.telescopes])
    n_tel = len(radii)

    # The sweeps run against a subsample; the winner is rescored against all of
    # it. Both the candidate pointings and the scored directions come out of the
    # region itself, so neither can wander off it.
    scored = (rng.choice(n_pixels, search_pixels, replace=False)
              if n_pixels > search_pixels else np.arange(n_pixels))
    aimed = (rng.choice(n_pixels, candidates, replace=False)
             if n_pixels > candidates else np.arange(n_pixels))

    sky = directions[scored]
    w = weights[scored] / weights[scored].sum()
    aims = directions[aimed]

    # Which scored directions each candidate pointing would cover, one matrix per
    # distinct camera. float32 so the marginal-gain step is a BLAS call rather
    # than a boolean promotion, which is most of the runtime.
    covers = {radius: ((aims @ sky.T) >= np.cos(radius)).astype(np.float32)
              for radius in np.unique(radii)}

    def climb(target, placed):
        """Sweep the telescopes until none of them wants to move."""
        placed = np.asarray(placed).copy()
        multiplicity = np.zeros(len(sky), dtype=np.float32)
        for i in range(n_tel):
            multiplicity += covers[radii[i]][placed[i]]

        for _ in range(sweeps):
            moved = 0
            for i in range(n_tel):
                without = multiplicity - covers[radii[i]][placed[i]]
                gain = (target * (_utility(without + 1, m_cut, eps)
                                  - _utility(without, m_cut, eps))).astype(np.float32)
                choice = int(np.argmax(covers[radii[i]] @ gain))
                if choice != placed[i]:
                    moved += 1
                    placed[i] = choice
                multiplicity = without + covers[radii[i]][placed[i]]
            if moved == 0:
                break

        value = float((target * _utility(multiplicity, m_cut, eps)).sum())
        return placed, multiplicity, value

    def build(target):
        """
        A first arrangement, telescopes laid down one at a time, widest first.

        Built against `m_cut` of one rather than the real one. With the real cut
        the first telescope gains nothing wherever it goes -- one camera never
        reaches two -- and the build has no gradient to follow. Covering the
        region first and letting the sweeps deepen it afterwards does.
        """
        placed = np.zeros(n_tel, dtype=int)
        multiplicity = np.zeros(len(sky), dtype=np.float32)
        for i in np.argsort(-radii):
            gain = (target * (_utility(multiplicity + 1, 1, eps)
                              - _utility(multiplicity, 1, eps))).astype(np.float32)
            placed[i] = int(np.argmax(covers[radii[i]] @ gain))
            multiplicity += covers[radii[i]][placed[i]]
        return placed

    def nearest(vectors):
        """The candidate pointings closest to an arrangement already made."""
        return np.argmax(aims @ np.asarray(vectors).T, axis=0)

    starts = []
    for count in warm_starts:
        if count > n_tel:
            continue
        point_subarrays(array, region, count, m_cut=m_cut)
        starts.append(nearest(array.pointing_vectors))

    available = reach(array, region_area)
    budget = available * len(sky)
    best = None
    for _ in range(rounds):
        target = target_field(w, budget, beta, m_cut)
        target = target / target.sum()

        best = None
        for start in [build(target)] + starts:
            trial = climb(target, start)
            if best is None or trial[2] > best[2]:
                best = trial

        # What the arrangement actually put on the region, which is what the
        # next target should be built from.
        landed = float(best[1].sum())
        starts = [best[0]] + starts[:2]
        if abs(landed - budget) <= 0.02 * max(budget, 1e-9):
            budget = landed
            break
        budget = 0.5 * (budget + landed)

    _point(array, aims[best[0]])
    return {
        "budget": budget / len(sky),
        "reach": available,
        "efficiency": budget / len(sky) / available,
        "objective": best[2],
        "pointings": [{"alt": alt, "az": az}
                      for alt, az in (as_altaz(v) for v in aims[best[0]])],
    }


def shaped_pointing_by_type(array, region, region_area=None, beta=1.0, m_cut=2,
                            **kwargs):
    """
    Shape each telescope type's own coverage, rather than the array's as one.

    `shaped_pointing` puts every telescope into one coordinate ascent with one
    shared objective, which lets a narrow-camera telescope and a wide one trade
    places freely -- and an SST and an MST are not fungible that way, any more
    than a sub-array is well served by mixing them. This solves the same
    objective once per telescope type instead, each with only its own cameras to
    place.

    `region_area` is passed unscaled to every type's solve. `reach` is camera
    solid angle divided by region area, so it is already per-type once each
    type's own camera radii go in; dividing the region between the types first
    would be answering a different question, since it is the same region every
    type is being asked to cover, not a slice of it.

    Parameters
    ----------
    array: `divtel.telescope.Array`
        left pointed at the combined answer
    region: `divtel.region.SkyRegion`
    region_area: `astropy.Quantity`, optional
    beta, m_cut: as `shaped_pointing`
    **kwargs:
        passed through to `shaped_pointing` for every type

    Returns
    -------
    dict
        ``types``, one entry per field-of-view radius holding that type's own
        `shaped_pointing` result; everything `describe` reports, measured on the
        whole array once every type is pointed.
    """
    if region_area is None:
        region_area = region.area

    types = {name: shaped_pointing(group, region, region_area, beta, m_cut,
                                   **kwargs)
             for name, group in array.group_by("fov_radius").items()}

    result = {"types": types}
    result.update(describe(array, region, m_cut))
    return result


def describe(array, region, m_cut=2):
    """
    What an arrangement delivers over a region.

    Parameters
    ----------
    array: `divtel.telescope.Array`
        already pointed, by anything -- this measures a divergent configuration,
        a set of sub-arrays or a shaped one alike
    region: `divtel.region.SkyRegion`
    m_cut: int

    Returns
    -------
    dict
        ``covered`` and ``covered_stereo``, the weight seen by one camera and by
        `m_cut`; ``blind``, the weight no camera sees; ``mean``, the
        weight-weighted mean multiplicity over the whole region; ``mean_covered``,
        the same average taken over the covered part alone; ``max``; and
        ``mismatch``.

        The two means differ only when coverage is incomplete, and then by a lot:
        an array covering 96 per cent of a region at depth 8.5 has a whole-region
        mean of 8.2, because the 4 per cent it abandoned enters as zero. Quote
        ``mean_covered`` beside ``covered_stereo`` and neither can flatter the
        other.

    Notes
    -----
    ``mismatch`` is the Kullback-Leibler divergence from the weight to the
    normalised multiplicity, over the part of the region that clears `m_cut`. It
    is zero when multiplicity is exactly proportional to weight there, and it is
    the number `shaped_pointing` is trying to make small. Read it beside
    ``covered_stereo`` and never alone: an array that abandons everything but the
    peak matches the shape of what it kept, perfectly, and is useless.

    It does not reach zero in practice. Multiplicity is a sum of camera discs, so
    as a field on the sky it cannot be sharper than one camera, and a region with
    structure finer than that -- which a two-detector arc has, across its width
    -- cannot be matched however the telescopes are placed. `camera_blur`
    measures that scale, so read the mismatch against it: the camera is the
    resolution limit of this whole idea.
    """
    weights = np.asarray(region.weights, dtype=float)
    weights = weights / weights.sum()

    multiplicity = region.multiplicity(array)
    stereo = multiplicity >= m_cut

    if stereo.any() and multiplicity[stereo].sum():
        p = weights[stereo] / weights[stereo].sum()
        q = multiplicity[stereo] / multiplicity[stereo].sum()
        mismatch = float((p * np.log(p / q)).sum())
    else:
        mismatch = float("nan")

    return {
        "covered": float(weights[multiplicity >= 1].sum()),
        "covered_stereo": float(weights[stereo].sum()),
        "blind": float(weights[multiplicity < 1].sum()),
        "mean": float((weights * multiplicity).sum()),
        "mean_covered": (float(np.average(multiplicity[stereo],
                                          weights=weights[stereo]))
                         if stereo.any() else 0.0),
        "max": int(multiplicity.max()),
        "mismatch": mismatch,
    }


def camera_blur(array, region, chunk=4096):
    """
    How far the weight map moves when one camera is smeared over it.

    A yardstick for the ``mismatch`` in `describe`, in the same units. It is not
    a bound and arrangements do beat it -- what it measures is the scale of the
    problem.

    Multiplicity is a sum of camera discs laid on the sky, so as a field it is
    the density of pointings smeared by a disc. Structure finer than a camera is
    therefore not reproducible, and a two-detector arc a few degrees wide against
    a several-degree camera has plenty of it. Blurring the map with one average
    camera and measuring how far the result has moved says how much of the
    mismatch is that limit rather than a failure of the pointing.

    A real array beats the single-blur number when it has a mix of cameras -- the
    narrow ones can go where the map is sharp -- so treat this as the scale at
    which matching the map stops being about where the telescopes point and
    starts being about how big their cameras are.

    Parameters
    ----------
    array: `divtel.telescope.Array`
        used for its camera sizes only; it is not re-pointed
    region: `divtel.region.SkyRegion`
    chunk: int
        directions smoothed at a time, so the pairwise comparison never has to
        exist all at once

    Returns
    -------
    float
        divergence from the map to its blurred self, comparable with
        `describe`'s ``mismatch``
    """
    directions = region.directions
    weights = np.asarray(region.weights, dtype=float)
    weights = weights / weights.sum()

    radii = np.array([t.fov_radius.to_value(u.rad) for t in array.telescopes])
    solid = 2 * np.pi * (1 - np.cos(radii))
    threshold = np.cos(float(np.average(radii, weights=solid)))

    blurred = np.empty(len(directions))
    for start in range(0, len(directions), chunk):
        block = directions[start:start + chunk]
        blurred[start:start + chunk] = ((block @ directions.T) >= threshold) @ weights

    blurred = blurred / blurred.sum()
    good = blurred > 0
    return float((weights[good] * np.log(weights[good] / blurred[good])).sum())


def multiplicity_by_probability(array, region, bins=5):
    """
    Mean multiplicity in each fifth of the weight, least likely first.

    The plain reading of whether the array is looking hardest where the source
    probably is. Directions are ordered by weight and cut into bins holding equal
    *total weight*, not equal area, so the last bin is the small bright core and
    the first is the wide faint skirt. A flat list is an array spread evenly over
    the region; a rising one is the point of `shaped_pointing`.

    Parameters
    ----------
    array: `divtel.telescope.Array`
        already pointed
    region: `divtel.region.SkyRegion`
    bins: int

    Returns
    -------
    list of float
        mean multiplicity per bin, in order of increasing weight
    """
    weights = np.asarray(region.weights, dtype=float)
    weights = weights / weights.sum()
    multiplicity = region.multiplicity(array)

    order = np.argsort(weights)
    cumulative = np.cumsum(weights[order])
    edges = np.linspace(0.0, 1.0, bins + 1)

    profile = []
    for low, high in zip(edges[:-1], edges[1:], strict=True):
        chosen = order[(cumulative > low) & (cumulative <= high)]
        profile.append(float(np.average(multiplicity[chosen],
                                        weights=weights[chosen]))
                       if len(chosen) else 0.0)
    return profile

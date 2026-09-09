"""Covering a weighted region: tiling, sub-arrays, and shaped pointing."""

from importlib.resources import files

import astropy.units as u
import numpy as np
import pytest

from divtel import strategy
from divtel.layout import load_array
from divtel.pointing import as_altaz
from divtel.region import SkyRegion

DATA = files("divtel") / "data"
SOUTH = DATA / "cta-south-paranal-alpha-prod6.ecsv"
NORTH = DATA / "cta-north-lapalma-alpha-prod6.ecsv"

ALT, AZ, M_CUT = 60 * u.deg, 180 * u.deg, 2


@pytest.fixture
def array():
    return load_array(SOUTH)


@pytest.fixture(scope="module")
def first_alert():
    """
    The first GW170817 localization, placed over CTAO-South and cut at the
    horizon. The showcase case: an arc no divergence can cover.
    """
    whole = SkyRegion.from_table(DATA / "gw170817" / "bayestar_hl_90.ecsv.gz")
    region, _ = whole.place(ALT, AZ).visible_part(0 * u.deg)
    return region


@pytest.fixture(scope="module")
def good_map():
    """The three-detector map: compact, and nearly within reach unaided."""
    whole = SkyRegion.from_table(DATA / "gw170817" / "bayestar_hlv_90.ecsv.gz")
    region, _ = whole.place(ALT, AZ).visible_part(0 * u.deg)
    return region


# -- reach and the target field -------------------------------------------

def test_reach_falls_as_the_region_grows(array, first_alert, good_map):
    assert (strategy.reach(array, good_map.area)
            > strategy.reach(array, first_alert.area))


def test_target_field_spends_the_whole_budget():
    weights = np.random.default_rng(0).random(200)
    field = strategy.target_field(weights, budget=1000.0, beta=1.0, m_cut=2)

    assert field.sum() == pytest.approx(1000.0)
    assert (field >= 0).all()


def test_a_tight_budget_serves_the_likeliest_and_abandons_the_rest():
    """
    Below the floor there is not enough to go round, and the strategy says which
    part of the region it is giving up on rather than covering everything once.
    """
    weights = np.arange(1.0, 101.0)
    field = strategy.target_field(weights, budget=40.0, beta=1.0, m_cut=2)

    served = field > 0
    assert served.sum() == 20  # 40 units of budget at a floor of 2
    # And it is the twenty most probable directions that got them.
    assert served[-20:].all()
    assert not served[:-20].any()


# -- tiling ----------------------------------------------------------------

def test_tiling_never_loses_ground(array, first_alert):
    tiles = strategy.tile_region(array, first_alert, pointings=8, m_cut=M_CUT)
    cumulative = [tile["cumulative"] for tile in tiles]

    assert cumulative == sorted(cumulative)
    assert all(tile["new"] >= 0 for tile in tiles)
    assert cumulative[-1] > 0.9


def test_the_first_tile_is_the_greediest(array, first_alert):
    tiles = strategy.tile_region(array, first_alert, pointings=6, m_cut=M_CUT)
    assert tiles[0]["new"] == max(tile["new"] for tile in tiles)


# -- sub-arrays ------------------------------------------------------------

def test_a_split_is_a_partition(array):
    groups = strategy.split_into_subarrays(array, 4)

    ids = [t.id for group in groups.values() for t in group.telescopes]
    assert sorted(ids) == sorted(t.id for t in array.telescopes)
    sizes = [len(group.telescopes) for group in groups.values()]
    assert max(sizes) - min(sizes) <= 1


def test_a_split_cannot_leave_a_group_empty(array):
    with pytest.raises(ValueError, match="without leaving some empty"):
        strategy.split_into_subarrays(array, len(array.telescopes) + 1)


def test_sub_arrays_share_their_telescopes_with_the_parent(array):
    """Pointing a group points the parent's telescopes, which is what lets
    multiplicity be measured on the whole array afterwards."""
    groups = strategy.split_into_subarrays(array, 3)
    groups["1"].divergent_pointing(0.0, 42 * u.deg, 77 * u.deg)

    pointed = [t for t in array.telescopes if t.id in
               {g.id for g in groups["1"].telescopes}]
    assert all(t.alt.to_value(u.deg) == pytest.approx(42) for t in pointed)


def test_sub_arrays_beat_the_best_divergence_on_the_first_alert(array, first_alert):
    """
    The finding the study turns on. No divergence covers more than about 70 per
    cent of this map stereoscopically; a split of the same array covers all of
    it, in one exposure, with every telescope on target.
    """
    result = strategy.point_subarrays(array, first_alert, 14, M_CUT)

    assert result["covered"] > 0.95
    assert result["on_target"] == len(array.telescopes)
    assert result["multiplicity"] >= M_CUT


def test_more_groups_reach_further_until_they_run_out(array, first_alert):
    covered = [strategy.point_subarrays(load_array(SOUTH), first_alert, count,
                                        M_CUT)["covered"]
               for count in (2, 4, 8, 14)]
    assert covered[-1] > covered[0]


def test_splitting_by_type_gives_each_its_own_count(good_map):
    array = load_array(SOUTH)
    result = strategy.point_subarrays_by_type(array, good_map, m_cut=M_CUT)

    assert set(result["types"]) == set(array.group_by("fov_radius"))
    assert result["covered"] > 0.9


def test_a_type_too_small_to_reconstruct_is_reported_as_such(good_map):
    """A group of one cannot see anything stereoscopically, and says so rather
    than being handed a sub-array count that cannot help it."""
    array = load_array(NORTH)
    single = array.group_by({"lone": [1], "rest": range(2, 14)})["lone"]

    result = strategy.point_subarrays_by_type(single, good_map, m_cut=M_CUT)
    assert set(result["types"].values()) == {None}


# -- unequal splits --------------------------------------------------------

def test_a_weighted_split_is_unequal_and_still_a_partition(array, first_alert):
    result = strategy.weighted_split(array, first_alert, 8, M_CUT)

    assert sum(result["sizes"]) == len(array.telescopes)
    assert min(result["sizes"]) >= M_CUT
    assert max(result["sizes"]) > min(result["sizes"])


def test_only_the_unequal_split_produces_a_gradient(array, first_alert):
    """
    Even sub-arrays are flat by construction: the multiplicity inside a group's
    footprint is the group's size, and every group is the same size. Sizing the
    groups from the probability is what bends the profile.
    """
    even_array = load_array(SOUTH)
    strategy.point_subarrays(even_array, first_alert, 8, M_CUT)
    flat = strategy.multiplicity_by_probability(even_array, first_alert)

    strategy.weighted_split(array, first_alert, 8, M_CUT)
    rising = strategy.multiplicity_by_probability(array, first_alert)

    assert flat[-1] / flat[0] < 1.3
    assert rising[-1] / rising[0] > 1.8


def test_a_weighted_split_needs_enough_telescopes(array, first_alert):
    with pytest.raises(ValueError, match="at least"):
        strategy.weighted_split(array, first_alert, 40, M_CUT)


# -- shaped pointing -------------------------------------------------------

def test_shaped_pointing_tracks_the_probability(array, first_alert):
    """
    The claim of the whole module: multiplicity proportional to probability. The
    least likely fifth of the map should be seen shallowly and the most likely
    fifth deeply, with everything in between in order.
    """
    strategy.shaped_pointing(array, first_alert, beta=1.0, m_cut=M_CUT)
    profile = strategy.multiplicity_by_probability(array, first_alert)

    assert profile == sorted(profile)
    assert profile[-1] / profile[0] > 2
    assert strategy.describe(array, first_alert, M_CUT)["covered_stereo"] > 0.95


def test_the_mismatch_is_smallest_at_beta_one(array, first_alert):
    """
    A check on the derivation, not a tuning result. Maximising `sum w log m`
    gives `m` proportional to `w`, so asking for `w**beta` should match the map
    best at beta = 1 and worse either side.
    """
    mismatch = {}
    for beta in (0.0, 1.0, 2.0):
        candidate = load_array(SOUTH)
        strategy.shaped_pointing(candidate, first_alert, beta=beta, m_cut=M_CUT)
        mismatch[beta] = strategy.describe(candidate, first_alert,
                                           M_CUT)["mismatch"]

    assert mismatch[1.0] < mismatch[0.0]
    assert mismatch[1.0] < mismatch[2.0]


def test_shaped_pointing_repeats_itself(first_alert):
    one, two = load_array(SOUTH), load_array(SOUTH)
    first = strategy.shaped_pointing(one, first_alert, m_cut=M_CUT, seed=3)
    second = strategy.shaped_pointing(two, first_alert, m_cut=M_CUT, seed=3)

    assert first["objective"] == pytest.approx(second["objective"])
    assert np.allclose(one.pointing_vectors, two.pointing_vectors)


def test_shaped_pointing_leaves_every_telescope_somewhere_useful(array, first_alert):
    result = strategy.shaped_pointing(array, first_alert, m_cut=M_CUT)

    assert len(result["pointings"]) == len(array.telescopes)
    assert 0 < result["efficiency"] <= 1
    # Every telescope is aimed inside the region, never at the sky beyond it.
    assert (first_alert.directions @ array.pointing_vectors.T).max(axis=0).min() > 0


def test_shaped_by_type_keeps_the_types_apart(good_map):
    array = load_array(SOUTH)
    result = strategy.shaped_pointing_by_type(array, good_map, m_cut=M_CUT,
                                              warm_starts=(4,), rounds=1)

    assert set(result["types"]) == set(array.group_by("fov_radius"))
    assert result["covered_stereo"] > 0.9


# -- measurement -----------------------------------------------------------

def test_describe_separates_the_two_means(array, first_alert):
    """
    ``mean`` averages over the whole region and ``mean_covered`` over the part
    that clears the cut, so they part company exactly when coverage does.
    """
    strategy.point_subarrays(array, first_alert, 8, M_CUT)
    described = strategy.describe(array, first_alert, M_CUT)

    assert described["mean_covered"] >= described["mean"]
    assert described["covered"] >= described["covered_stereo"]
    assert described["blind"] == pytest.approx(1 - described["covered"])


def test_a_parallel_array_matches_the_shape_of_what_it_kept(array, first_alert):
    """
    Why ``mismatch`` must be read beside ``covered_stereo``. A parallel array
    piles every telescope on the peak and abandons the rest of the map; what it
    kept is the shape of its own footprint, so it scores well on a metric it has
    no business winning.
    """
    array.divergent_pointing(0.0, *as_altaz(first_alert.centroid))
    parallel = strategy.describe(array, first_alert, M_CUT)

    shaped_array = load_array(SOUTH)
    strategy.shaped_pointing(shaped_array, first_alert, m_cut=M_CUT)
    shaped = strategy.describe(shaped_array, first_alert, M_CUT)

    assert parallel["covered_stereo"] < 0.4 < shaped["covered_stereo"]
    # And on mismatch alone it is competitive, which is the trap.
    assert parallel["mismatch"] < 2 * shaped["mismatch"]


def test_camera_blur_is_the_scale_of_the_problem(array, first_alert):
    """Multiplicity is a sum of camera discs, so it cannot be sharper than one
    camera. The blur says how much structure is out of reach by construction."""
    assert 0 < strategy.camera_blur(array, first_alert) < 1


def test_quintiles_hold_equal_probability_not_equal_area(array, first_alert):
    array.divergent_pointing(0.03, ALT, AZ)
    profile = strategy.multiplicity_by_probability(array, first_alert, bins=5)

    assert len(profile) == 5
    assert all(np.isfinite(profile))

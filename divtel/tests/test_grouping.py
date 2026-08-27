"""An array of several telescope types splits into sub-arrays."""

from importlib.resources import files

import astropy.units as u
import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from divtel.layout import load_array  # noqa: E402
from divtel.telescope import Array, Telescope  # noqa: E402
from divtel.visualization import display_groups  # noqa: E402

LA_PALMA = files("divtel") / "data" / "la_palma_4LST_15MST.ecsv"

CTAO_TYPES = {"LST": range(1, 5), "MST": range(5, 20)}


@pytest.fixture
def la_palma():
    return load_array(LA_PALMA)


@pytest.fixture
def pointed(la_palma):
    la_palma.divergent_pointing(0.02, 70 * u.deg, 180 * u.deg)
    return la_palma


def test_grouping_by_ctao_ids_splits_lsts_from_msts(la_palma):
    groups = la_palma.group_by(CTAO_TYPES)

    assert list(groups) == ["LST", "MST"]
    assert [tel.id for tel in groups["LST"].telescopes] == [1, 2, 3, 4]
    assert [tel.id for tel in groups["MST"].telescopes] == list(range(5, 20))


def test_each_group_holds_one_kind_of_telescope(la_palma):
    groups = la_palma.group_by(CTAO_TYPES)

    assert {tel.focal for tel in groups["LST"].telescopes} == {28 * u.m}
    assert {tel.focal for tel in groups["MST"].telescopes} == {16 * u.m}


def test_grouping_by_camera_radius_finds_the_same_split(la_palma):
    """Telescopes of a type share a camera, so the radius alone recovers them."""
    by_radius = la_palma.group_by("camera_radius")
    by_id = la_palma.group_by(CTAO_TYPES)

    sizes = sorted(len(group.telescopes) for group in by_radius.values())
    assert sizes == [4, 15]

    lst = next(g for g in by_radius.values() if len(g.telescopes) == 4)
    assert [tel.id for tel in lst.telescopes] == [tel.id for tel in by_id["LST"].telescopes]


def test_a_group_is_an_array_with_its_own_barycenter(la_palma):
    groups = la_palma.group_by(CTAO_TYPES)

    for group in groups.values():
        expected = u.Quantity([tel.position for tel in group.telescopes]).mean(axis=0)
        assert np.allclose(group.barycenter.to_value(u.m), expected.to_value(u.m))

    # and it is not the barycenter of the whole array
    assert not np.allclose(groups["LST"].barycenter.to_value(u.m),
                           la_palma.barycenter.to_value(u.m))


def test_a_group_keeps_the_pointing_of_the_array_it_came_from(pointed):
    """Array.__init__ resets the pointing state; group_by has to carry it over."""
    groups = pointed.group_by(CTAO_TYPES)

    for group in groups.values():
        assert group.div == pointed.div
        assert group.mean_pointing == pointed.mean_pointing
        group.hyper_fov()  # would raise if the group looked unpointed


def test_a_group_sees_less_sky_than_the_whole_array(pointed):
    groups = pointed.group_by(CTAO_TYPES)

    whole, _ = pointed.hyper_fov()
    lst, _ = groups["LST"].hyper_fov()
    assert lst < whole


def test_groups_share_telescopes_so_repointing_the_array_repoints_them(la_palma):
    groups = la_palma.group_by(CTAO_TYPES)
    la_palma.divergent_pointing(0.03, 60 * u.deg, 90 * u.deg)

    assert groups["LST"].telescopes[0] is la_palma.telescopes[0]
    assert groups["LST"].telescopes[0].alt == la_palma.telescopes[0].alt


def test_groups_need_not_cover_the_array(la_palma):
    groups = la_palma.group_by({"the four LSTs": [1, 2, 3, 4]})

    assert len(groups) == 1
    assert len(groups["the four LSTs"].telescopes) == 4


def test_an_unknown_telescope_id_is_rejected(la_palma):
    with pytest.raises(ValueError, match="which the array does not have"):
        la_palma.group_by({"LST": [1, 2, 3, 4], "nope": [99]})


def test_a_telescope_in_two_groups_is_rejected(la_palma):
    with pytest.raises(ValueError, match="is in both group"):
        la_palma.group_by({"LST": [1, 2, 3, 4], "also LST": [4, 5]})


def test_an_unknown_grouping_rule_is_rejected(la_palma):
    with pytest.raises(ValueError, match="dict of ids or 'camera_radius'"):
        la_palma.group_by("focal")


def test_ids_are_telescope_ids_not_positions():
    """A hand-built array numbers from a shared counter, so ids need not start at 1."""
    telescopes = [Telescope(x * u.m, 0 * u.m, 0 * u.m, 28 * u.m, 1 * u.m, tel_id=tel_id)
                  for x, tel_id in [(0, 41), (10, 42)]]
    groups = Array(telescopes).group_by({"pair": [41, 42]})

    assert [tel.id for tel in groups["pair"].telescopes] == [41, 42]
    with pytest.raises(ValueError, match="which the array does not have"):
        Array(telescopes).group_by({"pair": [1, 2]})


def test_display_groups_draws_one_colour_and_one_barycenter_per_group(pointed):
    ax = display_groups(pointed.group_by(CTAO_TYPES))

    assert [text.get_text() for text in ax.texts] == ["LST", "MST"]
    # per group: telescope scatter, telescope quiver, barycenter scatter, barycenter quiver
    assert len(ax.collections) == 4 * 2


@pytest.mark.parametrize("projection", ["xy", "xz", "yz"])
def test_display_groups_accepts_every_projection(pointed, projection):
    ax = display_groups(pointed.group_by(CTAO_TYPES), projection=projection)

    assert ax.get_xlabel() and ax.get_ylabel()


def test_display_groups_rejects_a_bad_projection(pointed):
    with pytest.raises(ValueError, match="projection should be"):
        display_groups(pointed.group_by(CTAO_TYPES), projection="ab")


def test_display_groups_rejects_nothing_to_draw():
    with pytest.raises(ValueError, match="no groups to display"):
        display_groups({})

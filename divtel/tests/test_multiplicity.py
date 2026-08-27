"""How much sky the array covers, and how well it covers it."""

from importlib.resources import files

import astropy.units as u
import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from divtel.layout import load_array  # noqa: E402
from divtel.telescope import Array, Telescope  # noqa: E402
from divtel.visualization import multiplicity_plot  # noqa: E402

LA_PALMA = files("divtel") / "data" / "la_palma_4LST_15MST.ecsv"


@pytest.fixture(autouse=True)
def a_figure_of_its_own():
    """The plots default to `plt.gca()`, which would otherwise be shared."""
    plt.close("all")
    yield
    plt.close("all")


@pytest.fixture
def square_array():
    """Four identical telescopes on a 100 m square."""
    return Array([
        Telescope(x * u.m, y * u.m, 0 * u.m, 20 * u.m, 1 * u.m)
        for x, y in [(100, 0), (-100, 0), (0, 100), (0, -100)]
    ])


@pytest.fixture
def la_palma():
    return load_array(LA_PALMA)


def test_the_profile_accounts_for_all_the_covered_sky(square_array):
    square_array.divergent_pointing(0.05, 70 * u.deg, 180 * u.deg)
    area, _ = square_array.hyper_fov()

    _, per_multiplicity = square_array.multiplicity_profile()

    assert per_multiplicity.sum().to_value(u.deg ** 2) == pytest.approx(
        area.to_value(u.deg ** 2)
    )


def test_a_parallel_array_is_seen_by_every_telescope_everywhere(square_array):
    """div=0 stacks every camera on the same disc, so multiplicity is exactly n."""
    square_array.divergent_pointing(0, 70 * u.deg, 180 * u.deg)

    multiplicity, _ = square_array.multiplicity_profile()
    mean, variance = square_array.multiplicity_moments()

    assert list(multiplicity) == [4]
    assert mean == pytest.approx(4)
    assert variance == pytest.approx(0)


def test_a_fully_divergent_array_is_seen_by_one_telescope_everywhere(square_array):
    square_array.divergent_pointing(1, 70 * u.deg, 180 * u.deg)

    multiplicity, _ = square_array.multiplicity_profile()
    mean, variance = square_array.multiplicity_moments()

    assert list(multiplicity) == [1]
    assert mean == pytest.approx(1)
    assert variance == pytest.approx(0)


@pytest.mark.parametrize("div", [0, 0.01, 0.02, 0.05, 0.1])
def test_mean_multiplicity_falls_as_divergence_grows(square_array, div):
    """The trade divergence makes: more sky, fewer telescopes on any of it."""
    means = []
    for value in [0, 0.01, 0.02, 0.05, 0.1]:
        square_array.divergent_pointing(value, 70 * u.deg, 180 * u.deg)
        means.append(square_array.multiplicity_moments()[0])

    assert means == sorted(means, reverse=True)


def test_area_and_multiplicity_pull_against_each_other(square_array):
    areas, means = [], []
    for div in [0, 0.02, 0.05, 0.1]:
        square_array.divergent_pointing(div, 70 * u.deg, 180 * u.deg)
        areas.append(square_array.hyper_fov()[0].to_value(u.deg ** 2))
        means.append(square_array.multiplicity_moments()[0])

    assert areas == sorted(areas)
    assert means == sorted(means, reverse=True)


def test_a_mixed_array_pointed_in_parallel_splits_by_camera(la_palma):
    """
    The 15 MSTs have the wider camera, the 4 LSTs the narrower one. Pointed in
    parallel, the sky the LSTs see is inside the sky the MSTs see, so there are
    two multiplicities: 19 over the LST disc and 15 over the rest.
    """
    la_palma.divergent_pointing(0, 70 * u.deg, 180 * u.deg)

    multiplicity, _ = la_palma.multiplicity_profile()

    assert list(multiplicity) == [15, 19]


def test_the_variance_is_zero_only_when_every_patch_is_seen_alike(la_palma):
    la_palma.divergent_pointing(0, 70 * u.deg, 180 * u.deg)

    _, variance = la_palma.multiplicity_moments()

    assert variance > 0


def test_weighting_is_by_area_not_by_patch(la_palma):
    """A sliver of overlap must not count for as much sky as the whole field."""
    la_palma.divergent_pointing(0, 70 * u.deg, 180 * u.deg)

    multiplicity, area = la_palma.multiplicity_profile()
    mean, _ = la_palma.multiplicity_moments()

    unweighted = np.mean(multiplicity)
    weighted = np.average(multiplicity, weights=area.to_value(u.deg ** 2))
    assert mean == pytest.approx(weighted)
    assert mean != pytest.approx(unweighted)


def test_patches_may_be_handed_in_rather_than_recomputed(la_palma):
    la_palma.divergent_pointing(0.02, 70 * u.deg, 180 * u.deg)
    _, patches = la_palma.hyper_fov()

    assert la_palma.multiplicity_moments(patches=patches) == \
        la_palma.multiplicity_moments()
    reused, fresh = (la_palma.multiplicity_profile(patches=patches),
                     la_palma.multiplicity_profile())
    assert list(reused[0]) == list(fresh[0])
    assert np.allclose(reused[1].value, fresh[1].value)


def test_multiplicity_plot_draws_one_bar_per_multiplicity(la_palma):
    la_palma.divergent_pointing(0.05, 70 * u.deg, 180 * u.deg)
    multiplicity, _ = la_palma.multiplicity_profile()

    ax = multiplicity_plot(la_palma)

    assert len(ax.patches) == len(multiplicity)
    assert [bar.get_x() + bar.get_width() / 2 for bar in ax.patches] == \
        pytest.approx(list(multiplicity))


def test_multiplicity_plot_bar_heights_are_the_profile(la_palma):
    la_palma.divergent_pointing(0.05, 70 * u.deg, 180 * u.deg)
    _, area = la_palma.multiplicity_profile()

    ax = multiplicity_plot(la_palma)

    assert [bar.get_height() for bar in ax.patches] == \
        pytest.approx(list(area.to_value(u.deg ** 2)))


def test_multiplicity_plot_fades_the_bars_below_the_cut(la_palma):
    la_palma.divergent_pointing(0.05, 70 * u.deg, 180 * u.deg)

    ax = multiplicity_plot(la_palma, m_cut=3)

    faded = [bar.get_alpha() for bar in ax.patches[:2]]
    kept = [bar.get_alpha() for bar in ax.patches[2:]]
    assert all(alpha is not None and alpha < 1 for alpha in faded)
    assert all(alpha is None for alpha in kept)


def test_multiplicity_plot_reports_the_area_above_the_cut(la_palma):
    la_palma.divergent_pointing(0.05, 70 * u.deg, 180 * u.deg)
    above_cut, _ = la_palma.hyper_fov(m_cut=3)

    ax = multiplicity_plot(la_palma, m_cut=3)

    assert f"{above_cut.value:.1f}" in ax.get_title()
    assert r"seen by $\geq$3" in ax.get_title()

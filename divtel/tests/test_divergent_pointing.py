import astropy.units as u
import numpy as np
import pytest

from divtel import pointing
from divtel.telescope import Array, Telescope


@pytest.fixture
def square_array():
    """Four telescopes on a 100 m square."""
    return Array(
        [
            Telescope(100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            Telescope(0 * u.m, 100 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            Telescope(-100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            Telescope(0 * u.m, -100 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
        ]
    )


@pytest.mark.parametrize("div", [1.5, 2.0, -0.2, -1.0, np.nan])
def test_divergent_pointing_rejects_div_out_of_range(square_array, div):
    """div outside [0, 1] is meaningless and must not be silently accepted."""
    with pytest.raises(ValueError, match="div must be between 0 and 1"):
        square_array.divergent_pointing(div, 70 * u.deg, 180 * u.deg)


@pytest.mark.parametrize("div", [0, 0.1, 0.5, 1.0])
def test_divergent_pointing_accepts_div_in_range(square_array, div):
    """The documented range stays usable, and produces real pointings."""
    square_array.divergent_pointing(div, 70 * u.deg, 180 * u.deg)
    alts = square_array.pointing_altaz[:, 0].to_value(u.deg)
    assert np.all(np.isfinite(alts))


def test_divergent_pointing_stays_above_the_horizon(square_array):
    """A pointing requested above the horizon must not end up below it."""
    square_array.divergent_pointing(0.5, 70 * u.deg, 180 * u.deg)
    alts = square_array.pointing_altaz[:, 0].to_value(u.deg)
    assert np.all(alts > 0), f"telescopes pointing underground: {alts}"


def test_divergent_pointing_zero_div_is_parallel(square_array):
    """div = 0 is the parallel case: every telescope on the mean direction."""
    square_array.divergent_pointing(0, 70 * u.deg, 180 * u.deg)
    altaz = square_array.pointing_altaz.to_value(u.deg)
    np.testing.assert_allclose(altaz[:, 0], 70)
    np.testing.assert_allclose(altaz[:, 1], 180)


@pytest.fixture
def hilly_pair():
    """
    Two telescopes on a slope, far enough apart in height that divergence
    pushes one of them below the horizon before it is clamped.
    """
    return Array([
        Telescope(100 * u.m, 0 * u.m, 50 * u.m, 20 * u.m, 1 * u.m),
        Telescope(-100 * u.m, 0 * u.m, -50 * u.m, 20 * u.m, 1 * u.m),
    ])


def test_divergent_pointing_clamps_a_telescope_to_the_horizon(hilly_pair):
    """
    Without clamping this telescope would end up at alt=-3.06 deg (checked
    against `pointing.tel_div_pointing` directly): a real mount cannot follow
    the ground, so divergent_pointing rails it at the horizon instead.
    """
    hilly_pair.divergent_pointing(0.3, 5 * u.deg, 180 * u.deg)

    alts = hilly_pair.pointing_altaz[:, 0].to_value(u.deg)
    assert alts[1] == pytest.approx(0.0)
    assert alts[0] > 0  # the other telescope was never in danger; untouched


def test_the_clamp_only_touches_altitude_not_azimuth(hilly_pair):
    """Azimuth keeps whatever the geometry gave it; only altitude rails."""
    g_point = pointing.pointG_position(hilly_pair.barycenter, 0.3, 5 * u.deg, 180 * u.deg)
    _, raw_az = pointing.tel_div_pointing(hilly_pair.telescopes[1].position, g_point)

    hilly_pair.divergent_pointing(0.3, 5 * u.deg, 180 * u.deg)

    assert hilly_pair.telescopes[1].az.to_value(u.deg) == pytest.approx(
        raw_az.to_value(u.deg)
    )


def test_no_telescope_ever_ends_up_below_the_horizon(hilly_pair):
    """The general guarantee, swept across the range that produced #40's bug."""
    for div in (0.1, 0.3, 0.5, 0.75, 1.0):
        for alt_mean in (-10, 0, 5, 20):
            hilly_pair.divergent_pointing(div, alt_mean * u.deg, 180 * u.deg)
            alts = hilly_pair.pointing_altaz[:, 0].to_value(u.deg)
            assert np.all(alts >= -1e-9), (div, alt_mean, alts)


def test_a_negative_alt_mean_is_also_clamped_in_the_parallel_case(square_array):
    """div=0 points every telescope straight at alt_mean; that rails too."""
    square_array.divergent_pointing(0, -20 * u.deg, 180 * u.deg)

    alts = square_array.pointing_altaz[:, 0].to_value(u.deg)
    np.testing.assert_allclose(alts, 0.0)


def test_mean_pointing_records_what_was_asked_not_what_was_achieved(hilly_pair):
    """
    Clamping is a per-telescope thing; the array still remembers the
    direction it was asked to diverge from, which export_cfg needs.
    """
    hilly_pair.divergent_pointing(0.3, 5 * u.deg, 180 * u.deg)

    alt_mean, az_mean = hilly_pair.mean_pointing
    assert alt_mean.to_value(u.deg) == pytest.approx(5.0)

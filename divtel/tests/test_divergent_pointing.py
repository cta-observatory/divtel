import astropy.units as u
import numpy as np
import pytest

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

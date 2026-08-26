"""The geometry API carries astropy units rather than bare metres."""

import astropy.units as u
import numpy as np
import pytest

from divtel import pointing
from divtel.telescope import Array, Telescope


@pytest.fixture
def array():
    return Array(
        [
            Telescope(100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            Telescope(0 * u.m, 100 * u.m, 5 * u.m, 20 * u.m, 1 * u.m),
            Telescope(-100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            Telescope(0 * u.m, -100 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
        ]
    )


def test_position_is_a_quantity():
    """`position` must keep its units, not build and then discard them."""
    telescope = Telescope(10 * u.m, 20 * u.m, 5 * u.m, 28 * u.m, 1 * u.m)
    position = telescope.position

    assert isinstance(position, u.Quantity)
    assert position.unit.is_equivalent(u.m)
    np.testing.assert_allclose(position.to_value(u.m), [10., 20., 5.])


def test_position_honours_the_input_unit():
    """A telescope built in km reports the same place as one built in m."""
    in_km = Telescope(0.01 * u.km, 0.02 * u.km, 0.005 * u.km, 28 * u.m, 1 * u.m)
    np.testing.assert_allclose(in_km.position.to_value(u.m), [10., 20., 5.])


def test_positions_array_and_barycenter_are_quantities(array):
    """The array-level aggregates must not strip the units either."""
    assert isinstance(array.positions_array, u.Quantity)
    assert array.positions_array.unit.is_equivalent(u.m)
    assert array.positions_array.shape == (4, 3)

    assert isinstance(array.barycenter, u.Quantity)
    np.testing.assert_allclose(array.barycenter.to_value(u.m), [0., 0., 1.25])


def test_point_to_object_accepts_a_quantity():
    """An object given in km must point the same way as one given in m."""
    telescope = Telescope(0 * u.m, 0 * u.m, 0 * u.m, 28 * u.m, 1 * u.m)

    telescope.point_to_object(np.array([1000., 1000., 1000.]))
    from_bare_array = telescope.pointing_vector

    telescope.point_to_object([1., 1., 1.] * u.km)
    np.testing.assert_allclose(telescope.pointing_vector, from_bare_array, atol=1e-12)


def test_pointing_helpers_return_quantities(array):
    """pointG_position and tel_div_pointing speak units too."""
    g_point = pointing.pointG_position(array.barycenter, 0.5, 70 * u.deg, 180 * u.deg)
    assert isinstance(g_point, u.Quantity)
    assert g_point.unit.is_equivalent(u.m)

    alt, az = pointing.tel_div_pointing(array.telescopes[0].position, g_point)
    assert alt.unit.is_equivalent(u.rad)
    assert az.unit.is_equivalent(u.rad)


def test_divergent_pointing_is_unit_agnostic():
    """Building the array in km must not change where it ends up pointing."""

    def build(unit, scale):
        return Array(
            [
                Telescope(x * scale * unit, y * scale * unit, 0 * unit,
                          20 * u.m, 1 * u.m)
                for x, y in [(100, 0), (0, 100), (-100, 0), (0, -100)]
            ]
        )

    in_m = build(u.m, 1.0)
    in_km = build(u.km, 0.001)
    for a in (in_m, in_km):
        a.divergent_pointing(0.5, 70 * u.deg, 180 * u.deg)

    np.testing.assert_allclose(in_m.pointing_altaz.to_value(u.deg),
                               in_km.pointing_altaz.to_value(u.deg),
                               atol=1e-9)

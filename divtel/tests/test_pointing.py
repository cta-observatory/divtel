"""Geometry helpers in `divtel.pointing`, tested directly."""

import astropy.units as u
import numpy as np
import pytest

from divtel import pointing

ALT_AZ_GRID = [
    (90, 0), (70, 180), (45, 90), (20, 270), (70, 30), (28, 0), (25, 0),
    (0, 45), (-30, 200), (89.999, 123),
]


@pytest.mark.parametrize("alt,az", ALT_AZ_GRID)
def test_local_frame_is_orthonormal_and_tangent(alt, az):
    increasing_az, increasing_alt = pointing.local_frame(alt * u.deg, az * u.deg)
    centre = pointing.alt_az_to_vector(alt * u.deg, az * u.deg)

    assert increasing_az @ increasing_az == pytest.approx(1.0)
    assert increasing_alt @ increasing_alt == pytest.approx(1.0)
    assert increasing_az @ increasing_alt == pytest.approx(0.0, abs=1e-9)
    assert centre @ increasing_az == pytest.approx(0.0, abs=1e-9)
    assert centre @ increasing_alt == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("alt,az", ALT_AZ_GRID)
def test_increasing_alt_moves_along_its_own_axis_only(alt, az):
    """A small step in altitude should land purely on the `increasing_alt` axis."""
    increasing_az, increasing_alt = pointing.local_frame(alt * u.deg, az * u.deg)
    before = pointing.alt_az_to_vector(alt * u.deg, az * u.deg)
    after = pointing.alt_az_to_vector((alt + 0.01) * u.deg, az * u.deg)
    step = after - before

    assert step @ increasing_alt > 0
    assert step @ increasing_az == pytest.approx(0.0, abs=1e-8)


@pytest.mark.parametrize("alt,az", ALT_AZ_GRID)
def test_increasing_az_moves_along_its_own_axis_only(alt, az):
    """A small step in azimuth should land purely on the `increasing_az` axis."""
    increasing_az, increasing_alt = pointing.local_frame(alt * u.deg, az * u.deg)
    before = pointing.alt_az_to_vector(alt * u.deg, az * u.deg)
    after = pointing.alt_az_to_vector(alt * u.deg, (az + 0.01) * u.deg)
    step = after - before

    assert step @ increasing_az > 0
    assert step @ increasing_alt == pytest.approx(0.0, abs=1e-8)


def test_orientation_is_continuous_through_the_old_seed_switch_region():
    """
    Regression test: an earlier implementation picked the tangent basis from a
    basis fixed to an external axis, seeded from one of two fixed vectors
    depending on which the pointing was closer to. That produced a real
    discontinuity -- crossing alt ~= 25.8 deg at az=0 flipped the whole frame
    by 90 degrees for an infinitesimal change in pointing. `local_frame` is
    built from the pointing's own alt/az instead, so it has no such seam.
    """
    altitudes = np.linspace(20, 35, 200)
    directions = []
    for alt in altitudes:
        increasing_az, increasing_alt = pointing.local_frame(alt * u.deg, 0 * u.deg)
        # angle of the "az" axis relative to a fixed reference, so a jump
        # would show up as a discontinuous step in this series
        directions.append(np.arctan2(increasing_az[1], increasing_az[0]))

    steps = np.abs(np.diff(np.unwrap(directions)))
    assert steps.max() < np.radians(1.0)


def test_local_frame_is_well_defined_exactly_at_the_zenith():
    """No division by cos(alt): the frame is still orthonormal at the pole."""
    increasing_az, increasing_alt = pointing.local_frame(90 * u.deg, 37 * u.deg)

    assert np.all(np.isfinite(increasing_az))
    assert np.all(np.isfinite(increasing_alt))
    assert increasing_az @ increasing_alt == pytest.approx(0.0, abs=1e-9)

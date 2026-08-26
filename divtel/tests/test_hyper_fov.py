import astropy.units as u
import numpy as np
import pytest

from divtel.telescope import Array, Telescope


@pytest.fixture
def hess_1():
    """Four telescopes on a 100 m square, all pointing at the zenith."""
    array = Array(
        [
            Telescope(100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            Telescope(0 * u.m, 100 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            Telescope(-100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            Telescope(0 * u.m, -100 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
        ]
    )
    for telescope in array.telescopes:
        telescope.point_to_altaz(90 * u.deg, 0 * u.deg)
    return array


def single_disc_area(array):
    radius = array.telescopes[0].fov_radius.to_value(u.deg)
    return np.pi * radius**2


def test_fov_radius_matches_fov_area():
    """The angular radius and the field of view area describe the same cone."""
    telescope = Telescope(0 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m)
    radius = telescope.fov_radius.to_value(u.rad)
    # `fov` uses the small-angle form, so agreement is close but not exact.
    assert np.pi * radius**2 == pytest.approx(telescope.fov.to_value(u.rad**2),
                                              rel=1e-2)


def test_parallel_pointing_is_one_camera_seen_by_all(hess_1):
    """With no divergence every disc coincides: one patch, full multiplicity."""
    hess_1.divergent_pointing(0, 70 * u.deg, 0 * u.deg)
    area, patches = hess_1.hyper_fov(m_cut=1)

    assert len(patches) == 1
    assert patches[0][1] == len(hess_1.telescopes)
    assert area.to_value(u.deg**2) == pytest.approx(single_disc_area(hess_1),
                                                    rel=1e-3)


def test_full_divergence_is_every_camera_seen_alone(hess_1):
    """Fully divergent, the discs are disjoint and nothing is seen twice."""
    hess_1.divergent_pointing(1, 70 * u.deg, 0 * u.deg)
    area, patches = hess_1.hyper_fov(m_cut=1)
    n = len(hess_1.telescopes)

    assert len(patches) == n
    assert {m for _, m in patches} == {1}
    assert area.to_value(u.deg**2) == pytest.approx(n * single_disc_area(hess_1),
                                                    rel=1e-3)
    # Nothing overlaps, so there is no stereoscopic coverage at all.
    assert hess_1.hyper_fov(m_cut=2)[0].to_value(u.deg**2) == pytest.approx(0)


def test_divergence_widens_coverage(hess_1):
    """Coverage grows with divergence, between the two limits above."""
    areas = []
    for div in (0.0, 0.01, 0.02, 0.05):
        hess_1.divergent_pointing(div, 70 * u.deg, 0 * u.deg)
        areas.append(hess_1.hyper_fov(m_cut=1)[0].to_value(u.deg**2))

    assert areas == sorted(areas)
    assert areas[0] == pytest.approx(single_disc_area(hess_1), rel=1e-3)
    assert areas[-1] <= 4 * single_disc_area(hess_1) * 1.001


def test_m_cut_selects_a_subset(hess_1):
    """A higher multiplicity cut can only ever count less sky."""
    hess_1.divergent_pointing(0.008, 70 * u.deg, 0 * u.deg)

    covered = hess_1.hyper_fov(m_cut=1)[0]
    stereo = hess_1.hyper_fov(m_cut=2)[0]
    quadruple = hess_1.hyper_fov(m_cut=4)[0]

    assert covered > stereo > quadruple > 0 * u.deg**2

    # The cut changes the area reported, never the patches returned.
    assert len(hess_1.hyper_fov(m_cut=1)[1]) == len(hess_1.hyper_fov(m_cut=4)[1])


def test_patches_partition_the_covered_area(hess_1):
    """The patches tile the union: no gaps, no double counting."""
    hess_1.divergent_pointing(0.01, 70 * u.deg, 0 * u.deg)
    area, patches = hess_1.hyper_fov(m_cut=1)

    assert sum(p.area for p, _ in patches) == pytest.approx(
        area.to_value(u.deg**2)
    )
    for patch, _ in patches:
        for other, _ in patches:
            if patch is not other:
                assert patch.intersection(other).area == pytest.approx(0, abs=1e-9)


def test_azimuth_wrap_is_unwrapped(hess_1):
    """An array straddling due north is not torn in half by the 0/360 seam."""
    for telescope in hess_1.telescopes:
        telescope.point_to_altaz(70 * u.deg, 0 * u.deg)
    # Nudge two telescopes either side of the seam.
    hess_1.telescopes[0].point_to_altaz(70 * u.deg, 359 * u.deg)
    hess_1.telescopes[1].point_to_altaz(70 * u.deg, 1 * u.deg)

    area, patches = hess_1.hyper_fov(m_cut=1)
    # Those two discs are 2 deg apart and 5.7 deg across, so they must overlap.
    assert max(m for _, m in patches) >= 2
    assert area.to_value(u.deg**2) < 4 * single_disc_area(hess_1)

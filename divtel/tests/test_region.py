"""A weighted region of sky: how it is placed, cut, and measured."""

from importlib.resources import files

import astropy.units as u
import numpy as np
import pytest

from divtel.layout import load_array
from divtel.pointing import as_altaz
from divtel.region import SkyRegion

GW170817 = files("divtel") / "data" / "gw170817"
SOUTH = files("divtel") / "data" / "cta-south-paranal-alpha-prod6.ecsv"


def cap(radius, n=400, seed=0):
    """A disc of `radius` about the x axis, sampled evenly, uniformly weighted."""
    rng = np.random.default_rng(seed)
    cosines = 1 - rng.random(n) * (1 - np.cos(radius.to_value(u.rad)))
    angles = rng.random(n) * 2 * np.pi
    sines = np.sqrt(1 - cosines**2)
    directions = np.column_stack([cosines,
                                  sines * np.cos(angles),
                                  sines * np.sin(angles)])
    # Put one direction exactly on the rim, so the containment radius has
    # something to find rather than depending on how the sample fell.
    directions[0] = [np.cos(radius.to_value(u.rad)),
                     np.sin(radius.to_value(u.rad)), 0.0]
    return SkyRegion(directions, np.ones(n), pixel_area=1e-3 * u.deg**2)


def test_length_must_match():
    with pytest.raises(ValueError, match="against"):
        SkyRegion(np.zeros((3, 3)), np.ones(2))


def test_placement_is_rigid():
    """Placing a region moves it without changing any angle inside it."""
    region = cap(8 * u.deg)
    before = region.directions @ region.directions.T

    placed = region.place(35 * u.deg, 210 * u.deg)

    assert np.allclose(placed.directions @ placed.directions.T, before, atol=1e-12)
    assert np.allclose(np.linalg.norm(placed.directions, axis=1), 1.0)


def test_placement_puts_the_centroid_where_asked():
    alt, az = 35 * u.deg, 210 * u.deg
    placed = cap(8 * u.deg).place(alt, az)

    got_alt, got_az = as_altaz(placed.centroid)
    assert got_alt.to_value(u.deg) == pytest.approx(alt.to_value(u.deg), abs=1e-6)
    # divtel reports azimuth in (-180, 180], so 210 comes back as -150.
    assert got_az.to_value(u.deg) % 360 == pytest.approx(
        az.to_value(u.deg) % 360, abs=1e-6)


def test_a_shared_rotation_keeps_nested_bands_nested():
    """Two bands of one map move together when they share a rotation."""
    outer = cap(8 * u.deg)
    inner = outer[np.arange(0, len(outer), 3)]

    rotation = outer.placement_rotation(50 * u.deg, 100 * u.deg)
    moved_outer, moved_inner = outer.rotate(rotation), inner.rotate(rotation)

    # The inner band's directions are a subset of the outer one's, and they are
    # still exactly where the outer band left them.
    assert np.allclose(moved_inner.directions,
                       moved_outer.directions[np.arange(0, len(outer), 3)])


def test_containment_radius_of_a_cap_is_its_radius():
    region = cap(6 * u.deg)
    assert region.containment_radius().to_value(u.deg) == pytest.approx(6.0, abs=0.01)


def test_max_separation_of_a_cap_is_its_diameter():
    region = cap(6 * u.deg, n=3000)
    assert region.max_separation.to_value(u.deg) == pytest.approx(12.0, abs=0.5)


def test_equivalent_radius_needs_a_pixel_area():
    region = SkyRegion(cap(6 * u.deg).directions, np.ones(400))
    with pytest.raises(ValueError, match="no pixel area"):
        _ = region.area


def test_visible_part_keeps_the_weights_it_did_not_cut():
    """
    The survivors' weights are left alone, so coverage measured after the cut is
    still a fraction of the whole region rather than of what is left of it.
    """
    region = cap(40 * u.deg).place(20 * u.deg, 180 * u.deg)
    above, dropped = region.visible_part(0 * u.deg)

    assert 0 < dropped < 1
    assert (above.directions[:, 2] >= 0).all()
    assert above.weights.sum() == pytest.approx((1 - dropped) * region.weights.sum())
    assert len(above) + round(dropped * len(region)) == pytest.approx(len(region), abs=1)


def test_multiplicity_of_a_parallel_array_is_all_or_nothing():
    """
    Pointed normally, a direction is seen by whichever cameras reach it -- so
    inside the narrowest camera every telescope sees it, and past the widest
    none does.
    """
    array = load_array(SOUTH)
    array.divergent_pointing(0.0, 60 * u.deg, 180 * u.deg)
    radii = sorted(t.fov_radius for t in array.telescopes)

    inside = cap(0.5 * radii[0]).place(60 * u.deg, 180 * u.deg)
    outside = cap(3 * radii[-1]).place(60 * u.deg, 180 * u.deg)

    assert (inside.multiplicity(array) == len(array.telescopes)).all()
    assert outside.multiplicity(array).min() == 0


def test_covered_is_a_fraction_of_the_whole_region():
    array = load_array(SOUTH)
    array.divergent_pointing(0.0, 60 * u.deg, 180 * u.deg)
    region = cap(20 * u.deg).place(60 * u.deg, 180 * u.deg)

    assert 0 < region.covered(array, m_cut=2) < 1
    # More telescopes are never needed to see less.
    assert region.covered(array, m_cut=1) >= region.covered(array, m_cut=2)


def test_bundled_regions_reproduce_the_published_areas():
    """
    The 90% areas LIGO/Virgo published for the three GW170817 alert maps, as
    circulated in GCN 21509, 21513 and 21527.
    """
    published = {"bayestar_hl": 190, "bayestar_hlv": 31,
                 "lalinference_prelim": 33.6}

    for name, area in published.items():
        region = SkyRegion.from_table(GW170817 / f"{name}_90.ecsv.gz")
        assert region.area.to_value(u.deg**2) == pytest.approx(area, rel=0.02)
        assert region.weights.sum() == pytest.approx(0.9, abs=0.001)


def test_the_first_alert_map_is_an_arc_and_not_a_disc():
    """
    Area is not the same as the shape of the map. The first
    GW170817 map covers 187 deg2, which sounds like a patch 15 degrees across.
    It is 124 degrees long.
    """
    region = SkyRegion.from_table(GW170817 / "bayestar_hl_90.ecsv.gz")
    extent = region.describe_extent()

    assert extent["equivalent_radius"].to_value(u.deg) == pytest.approx(7.7, abs=0.1)
    assert extent["containment_radius"].to_value(u.deg) == pytest.approx(62.1, abs=0.2)
    assert extent["max_separation"].to_value(u.deg) == pytest.approx(124.3, abs=0.2)
    assert extent["elongation"] > 5


def test_region_survives_a_write_and_a_read(tmp_path):
    region = SkyRegion.from_table(GW170817 / "bayestar_hlv_90.ecsv.gz")
    path = tmp_path / "round_trip.ecsv.gz"
    region.to_table(path)

    back = SkyRegion.from_table(path)
    assert len(back) == len(region)
    assert back.meta["gcn"] == region.meta["gcn"]
    assert back.pixel_area == region.pixel_area
    # Written to five decimals of a degree, which is a tenth of an arcsecond.
    assert np.allclose(back.directions, region.directions, atol=1e-6)

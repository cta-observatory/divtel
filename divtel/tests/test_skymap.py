"""Reading a published HEALPix localization. Needs the `skymap` extra."""

from importlib.resources import files

import astropy.units as u
import numpy as np
import pytest

pytest.importorskip("astropy_healpix",
                    reason="divtel.skymap needs the divtel[skymap] extra")

from divtel import skymap  # noqa: E402
from divtel.layout import load_array  # noqa: E402
from divtel.pointing import alt_az_to_vector  # noqa: E402

DATA = files("divtel") / "data"
# The FITS maps are downloaded by data/gw170817/make_regions.py and not
# committed, so anything reading one is skipped where they are absent.
MAPS = DATA / "gw170817" / "_maps"

needs_maps = pytest.mark.skipif(
    not (MAPS / "bayestar_no_virgo.fits.gz").exists(),
    reason="run divtel/data/gw170817/make_regions.py to fetch the sky maps")


@pytest.fixture(scope="module")
def first_alert():
    path = skymap.fetch("bayestar_hl", MAPS, download=False)
    return skymap.load(path, skymap.SKYMAPS["bayestar_hl"], max_nside=128)


def test_every_published_map_carries_its_provenance():
    """A follow-up is planned against the best localization available when the
    telescope could move, so when each map arrived is part of the record."""
    for info in skymap.SKYMAPS.values():
        assert info.gcn and info.t_available and info.detectors
        assert info.url.startswith("https://")


def test_an_unknown_map_says_which_ones_it_knows(tmp_path):
    with pytest.raises(KeyError, match="bayestar_hl"):
        skymap.fetch("no_such_map", tmp_path)


def test_a_missing_map_is_an_error_when_downloading_is_off(tmp_path):
    with pytest.raises(FileNotFoundError, match="not cached"):
        skymap.fetch("bayestar_hl", tmp_path, download=False)


def test_credible_levels_order_the_pixels_by_probability():
    prob = np.array([0.5, 0.1, 0.3, 0.1])
    levels = skymap.credible_levels(prob)

    # The most probable pixel joins first, at its own probability.
    assert levels[0] == pytest.approx(0.5)
    assert levels[2] == pytest.approx(0.8)
    assert levels.max() == pytest.approx(1.0)


@needs_maps
def test_degrading_conserves_probability_exactly(first_alert):
    """The four children of a NESTED pixel are adjacent in the array, so
    coarsening is a reshape and a sum -- and a sum, not an average."""
    coarse = first_alert.degrade(64)

    assert coarse.nside == 64
    assert coarse.prob.sum() == pytest.approx(first_alert.prob.sum())
    # The area is quantised to the coarse pixel and so does move. It moves most
    # for a thin region, which has a great deal of edge per unit of area: this
    # one is an arc a few degrees wide against pixels nearly a degree across.
    assert skymap.credible_area(coarse, 0.9).to_value(u.deg**2) == pytest.approx(
        skymap.credible_area(first_alert, 0.9).to_value(u.deg**2), rel=0.15)


@needs_maps
def test_degrading_only_ever_coarsens(first_alert):
    with pytest.raises(ValueError, match="only makes a map coarser"):
        first_alert.degrade(512)
    with pytest.raises(ValueError, match="power of two"):
        first_alert.degrade(100)


@needs_maps
def test_the_two_detector_map_has_two_lobes(first_alert):
    """
    With two detectors the timing leaves an ambiguity that puts a second lobe
    somewhere else entirely, and no array covers both at once. An observer picks
    a lobe, so a study has to be able to.
    """
    pieces = skymap.split_regions(first_alert, 0.9)

    assert len(pieces) >= 2
    shares = [float(first_alert.prob[piece].sum()) for piece in pieces]
    assert shares[0] == pytest.approx(0.61, abs=0.02)
    assert shares[1] == pytest.approx(0.29, abs=0.02)
    # The pieces are exactly the credible region, cut up.
    union = np.logical_or.reduce(pieces)
    assert (union == skymap.credible_mask(first_alert, 0.9)).all()


@needs_maps
def test_a_region_carries_the_map_it_came_from(first_alert):
    region = skymap.region(first_alert, 0.9)

    assert region.meta["gcn"] == "GCN 21509"
    assert region.meta["level"] == 0.9
    assert region.weights.sum() == pytest.approx(0.9, abs=0.001)
    assert region.area.to_value(u.deg**2) == pytest.approx(190, rel=0.03)


@needs_maps
def test_one_lobe_can_be_carried_on_alone(first_alert):
    biggest = skymap.split_regions(first_alert, 0.9)[0]
    region = skymap.region(first_alert, 0.9, mask=biggest)

    assert region.weights.sum() == pytest.approx(0.61, abs=0.02)
    assert len(region) < len(skymap.region(first_alert, 0.9))


def test_pixel_counting_agrees_with_the_polygon_arithmetic():
    """
    Cross-check on the whole coverage story. `SkyRegion.multiplicity` tests each
    direction against each camera on the sphere; `Array.hyper_fov` cuts the same
    footprints into polygons in an equal-area projection. They share no geometry,
    so agreement is real evidence rather than one bug reproducing itself.
    """
    from astropy_healpix import HEALPix, nside_to_pixel_area

    from divtel.region import SkyRegion

    array = load_array(DATA / "cta-north-lapalma-alpha-prod6.ecsv")
    array.divergent_pointing(0.04, 70 * u.deg, 180 * u.deg)

    # A uniform grid read straight as ground-frame directions: no sky and no
    # time, just every direction over the array, which is all an area needs.
    healpix = HEALPix(nside=256)
    lon, lat = healpix.healpix_to_lonlat(np.arange(healpix.npix))
    grid = SkyRegion(np.column_stack(alt_az_to_vector(lat.to(u.rad), lon.to(u.rad))),
                     np.ones(healpix.npix),
                     pixel_area=nside_to_pixel_area(256).to(u.deg**2))

    multiplicity = grid.multiplicity(array)
    pixel_area = grid.pixel_area.to_value(u.deg**2)

    for m_cut in (1, 2):
        counted = (multiplicity >= m_cut).sum() * pixel_area
        polygons = array.hyper_fov(min_telescopes=m_cut)[0].to_value(u.deg**2)
        assert counted == pytest.approx(polygons, rel=0.01)

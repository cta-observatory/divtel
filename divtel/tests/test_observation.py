"""
Pointing at a real source at a real time.

Everything here runs offline: no test resolves a source name or looks a site up
in astropy's registry, since both go out to the network.
"""

from importlib.resources import files

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

from divtel.layout import load_array
from divtel.observation import SITES, Observation, pointing_coord

CRAB = SkyCoord(ra=83.633 * u.deg, dec=22.015 * u.deg)
NIGHT = "2026-03-01T23:00:00"


@pytest.fixture
def observation():
    return Observation(site="north", time=NIGHT)


@pytest.fixture
def array():
    return load_array(files("divtel") / "data" / "la_palma_4LST_15MST.ecsv")


def test_a_source_pointed_at_is_the_source_read_back(observation, array):
    """
    The one that matters: divtel's azimuth and astropy's must agree, or an array
    aimed at a source would silently watch somewhere else. Point every telescope
    at the Crab, then ask the telescopes where they are looking.
    """
    alt, az = observation.altaz_of(CRAB)
    array.divergent_pointing(0, alt, az)

    pointings = pointing_coord(array, observation)

    assert pointings.separation(CRAB).max() < 1 * u.arcsec


def test_a_divergent_array_spreads_around_the_source(observation, array):
    alt, az = observation.altaz_of(CRAB)
    array.divergent_pointing(0.02, alt, az)

    separations = pointing_coord(array, observation).separation(CRAB)

    assert separations.min() > 0 * u.deg
    assert separations.max() < 10 * u.deg


def test_pointings_can_be_left_in_the_horizontal_frame(observation, array):
    alt, az = observation.altaz_of(CRAB)
    array.divergent_pointing(0.02, alt, az)

    horizontal = pointing_coord(array, observation, icrs=False)

    assert np.allclose(horizontal.alt.to_value(u.deg),
                       array.pointing_altaz[:, 0].to_value(u.deg))


def test_the_pole_sits_due_north_at_the_latitude_of_the_site(observation):
    """A physical check on the site: no trigonometry of ours is involved."""
    pole = SkyCoord(ra=0 * u.deg, dec=90 * u.deg)

    alt, az = observation.altaz_of(pole)

    latitude = observation.location.geodetic.lat
    assert abs(alt - latitude) < 0.2 * u.deg
    assert min(az.to_value(u.deg), 360 - az.to_value(u.deg)) < 0.5


def test_a_source_moves_across_the_sky_as_the_night_goes_on(observation):
    alt_now, az_now = observation.altaz_of(CRAB)
    alt_later, az_later = observation.after(3 * u.hour).altaz_of(CRAB)

    assert abs(alt_later - alt_now) > 1 * u.deg
    assert abs(az_later - az_now) > 1 * u.deg


def test_the_sky_the_array_sees_follows_the_source(observation, array):
    """Tracking a setting source costs coverage, because the array tips over."""
    areas = []
    for hours in (0, 3, 6):
        moment = observation.after(hours * u.hour)
        alt, az = moment.altaz_of(CRAB)
        array.divergent_pointing(0.02, alt, az)
        areas.append(array.hyper_fov()[0])

    assert len({area.value for area in areas}) == 3


def test_at_and_after_leave_the_original_alone(observation):
    later = observation.after(2 * u.hour)
    elsewhen = observation.at("2026-06-01T23:00:00")

    assert observation.time == Time(NIGHT, scale="utc")
    assert (later.time - observation.time).to_value(u.hour) == pytest.approx(2)
    assert elsewhen.time != observation.time
    assert later.location == observation.location


def test_a_negative_delta_goes_backwards(observation):
    assert observation.after(-30 * u.min).time.isot == "2026-03-01T22:30:00.000"


def test_both_ctao_sites_are_where_they_should_be():
    """Catches a transposed or sign-flipped coordinate, which would be silent."""
    north = SITES["north"].geodetic
    south = SITES["south"].geodetic

    assert 0 * u.deg < north.lat < 45 * u.deg      # northern hemisphere
    assert -30 * u.deg < north.lon < 0 * u.deg     # Atlantic, west of Greenwich
    assert -45 * u.deg < south.lat < 0 * u.deg     # southern hemisphere
    assert -90 * u.deg < south.lon < -45 * u.deg   # Chile
    for site in (north, south):
        assert 1500 * u.m < site.height < 3000 * u.m


def test_an_earth_location_is_taken_as_it_is():
    somewhere = EarthLocation.from_geodetic(
        lon=10 * u.deg, lat=45 * u.deg, height=100 * u.m
    )

    assert Observation(site=somewhere, time=NIGHT).location == somewhere


def test_site_names_are_not_case_sensitive():
    assert Observation(site="North", time=NIGHT).location == SITES["north"]


def test_an_unknown_site_says_so_without_going_to_the_network():
    with pytest.raises(ValueError, match="unknown site"):
        Observation(site="Mount Doom", time=NIGHT)


def test_a_site_that_is_not_a_name_or_a_location_is_rejected():
    with pytest.raises(ValueError, match="must be a name or an EarthLocation"):
        Observation(site=42, time=NIGHT)


def test_the_time_defaults_to_now():
    assert abs(Observation(site="north").time - Time.now()) < 10 * u.s


def test_the_sun_is_down_at_night_and_up_at_noon(observation):
    """get_moon is gone from astropy; these go through get_body."""
    assert observation.sun.alt < 0 * u.deg
    assert observation.at("2026-03-01T13:00:00").sun.alt > 0 * u.deg


def test_the_moon_is_somewhere(observation):
    assert -90 * u.deg <= observation.moon.alt <= 90 * u.deg


def test_the_frame_carries_the_place_and_the_time(observation):
    frame = observation.altaz

    assert frame.obstime == observation.time
    assert frame.location == observation.location

"""
When and where the array is observing from.

The geometry in `divtel.telescope` knows nothing about the sky. It works in a
frame fixed to the ground -- x north, y west, z up -- and an array pointed at
alt 70, az 180 stays pointed there whatever the hour. That is the right model
for the geometry, and it is not enough to plan an observation: a source rises
and sets, so the direction to point at depends on where you are and what time
it is.

This module supplies that missing half, and stays optional. Nothing in
`divtel.telescope` imports it, and an `Observation` never touches an `Array` --
it converts a sky position into the alt/az pair `Array.divergent_pointing`
already takes:

    obs = Observation(site="north", time="2026-03-01T23:00:00")
    alt, az = obs.altaz_of(SkyCoord(ra=83.633 * u.deg, dec=22.015 * u.deg))
    array.divergent_pointing(0.02, alt, az)

divtel's azimuth is astropy's: measured from north through east. No conversion
happens anywhere in here, and a test pins that down by pointing an array at a
source and reading its right ascension and declination back off the telescopes.
"""

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
from astropy.time import Time

__all__ = ["Observation", "SITES", "pointing_coord"]


# The CTAO array-centre reference positions, as used to configure the
# production air-shower simulations. They are built in rather than looked up so
# that the common case needs no network, and because astropy's site registry
# has no CTAO-South entry -- only Paranal, the VLT site a couple of kilometres
# away. For work that turns on the exact position of a particular telescope,
# pass your own EarthLocation.
SITES = {
    "north": EarthLocation.from_geodetic(
        lon=-17.8917 * u.deg, lat=28.7619 * u.deg, height=2158 * u.m
    ),
    "south": EarthLocation.from_geodetic(
        lon=-70.3167 * u.deg, lat=-24.6833 * u.deg, height=2147 * u.m
    ),
}


def _resolve_site(site):
    """An `EarthLocation` from a built-in name, an astropy site name, or itself."""
    if isinstance(site, EarthLocation):
        return site

    if not isinstance(site, str):
        raise ValueError(
            f"site must be a name or an EarthLocation, got {type(site).__name__}"
        )

    if site.lower() in SITES:
        return SITES[site.lower()]

    try:
        return EarthLocation.of_site(site)
    except Exception as error:
        raise ValueError(
            f"unknown site {site!r}; use one of {', '.join(sorted(SITES))}, a name "
            "astropy's site registry knows, or an EarthLocation. Note that looking a "
            f"name up in that registry needs network access ({error})"
        ) from error


class Observation:
    """
    A place and a time to observe from.

    Parameters
    ----------
    site: str or `astropy.coordinates.EarthLocation`
        ``"north"`` or ``"south"`` for the CTAO sites; any other name is looked
        up in astropy's site registry, which needs network access the first time
        it is used. An `EarthLocation` is taken as it is.
    time: str, `astropy.time.Time`, optional
        when the observation happens, UTC. Defaults to now.

    Attributes
    ----------
    location: `astropy.coordinates.EarthLocation`
    time: `astropy.time.Time`

    Examples
    --------
    >>> obs = Observation(site="north", time="2026-03-01T23:00:00")
    >>> crab = SkyCoord(ra=83.633 * u.deg, dec=22.015 * u.deg)
    >>> alt, az = obs.altaz_of(crab)
    >>> f"{alt:.2f}, {az:.2f}"
    '50.93 deg, 270.16 deg'
    """

    def __init__(self, site="north", time=None):
        self.location = _resolve_site(site)
        self.time = Time.now() if time is None else Time(time, scale="utc")

    def __repr__(self):
        lon, lat, height = self.location.geodetic
        return (f"Observation(lon={lon.deg:.4f} deg, lat={lat.deg:.4f} deg, "
                f"height={height.to_value(u.m):.0f} m, time={self.time.isot})")

    @property
    def altaz(self):
        """
        The horizontal frame at this place and time.

        Returns
        -------
        `astropy.coordinates.AltAz`
        """
        return AltAz(obstime=self.time, location=self.location)

    def at(self, time):
        """
        The same place at another time.

        Returns a new `Observation` rather than changing this one, so an array
        pointed from one cannot quietly fall out of step with it.

        Parameters
        ----------
        time: str or `astropy.time.Time`

        Returns
        -------
        `Observation`
        """
        return Observation(site=self.location, time=time)

    def after(self, delta):
        """
        The same place, later (or, for a negative delta, earlier).

        Parameters
        ----------
        delta: `astropy.Quantity`
            elapsed time, e.g. ``2 * u.hour``

        Returns
        -------
        `Observation`

        Examples
        --------
        >>> obs.after(-30 * u.min).time.isot
        '2026-03-01T22:30:00.000'
        """
        return self.at(self.time + delta)

    def altaz_of(self, target):
        """
        Where to point to see a target, in the alt/az `Array` takes.

        Parameters
        ----------
        target: `astropy.coordinates.SkyCoord` or str
            a sky position, or a name to resolve -- resolving a name queries a
            name server, so it needs network access.

        Returns
        -------
        (alt, az): tuple of `astropy.Quantity`, in degrees

        Notes
        -----
        The azimuth is returned as astropy measures it, from north through east,
        which is already what `Array.divergent_pointing` expects.

        Examples
        --------
        >>> alt, az = obs.altaz_of(SkyCoord(ra=83.633 * u.deg, dec=22.015 * u.deg))
        >>> array.divergent_pointing(0.02, alt, az)
        """
        if isinstance(target, str):
            target = SkyCoord.from_name(target)

        horizontal = target.transform_to(self.altaz)
        return horizontal.alt.to(u.deg), horizontal.az.to(u.deg)

    def body(self, name):
        """
        Where a solar system body is, in this frame.

        Parameters
        ----------
        name: str
            as `astropy.coordinates.get_body` takes it, e.g. ``"moon"``

        Returns
        -------
        `astropy.coordinates.SkyCoord`
            in the horizontal frame of this observation
        """
        return get_body(name, self.time, location=self.location).transform_to(self.altaz)

    @property
    def sun(self):
        """
        The Sun, in this frame -- negative altitude means night.

        Returns
        -------
        `astropy.coordinates.SkyCoord`
        """
        return self.body("sun")

    @property
    def moon(self):
        """
        The Moon, in this frame.

        Returns
        -------
        `astropy.coordinates.SkyCoord`
        """
        return self.body("moon")


def pointing_coord(array, observation, icrs=True):
    """
    Where each telescope of an array is pointing, on the sky.

    The inverse of the `Observation.altaz_of` trip: alt/az in, sky position out.
    With divergent pointing the telescopes look in different directions, so this
    is one coordinate per telescope, and it is what says which part of the sky
    each one is actually watching.

    Parameters
    ----------
    array: `Array`
    observation: `Observation`
    icrs: bool
        return right ascension and declination; otherwise leave the coordinates
        in the horizontal frame

    Returns
    -------
    `astropy.coordinates.SkyCoord`
        one entry per telescope, in the order of ``array.telescopes``

    Examples
    --------
    >>> alt, az = obs.altaz_of(crab)
    >>> array.divergent_pointing(0.02, alt, az)
    >>> pointing_coord(array, obs).separation(crab).to(u.deg).max()
    <Angle 4.3121837 deg>
    """
    altaz = array.pointing_altaz

    coord = SkyCoord(alt=altaz[:, 0], az=altaz[:, 1], frame=observation.altaz)

    return coord.icrs if icrs else coord

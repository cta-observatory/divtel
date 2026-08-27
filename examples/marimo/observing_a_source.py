import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")


@app.cell(hide_code=True)
async def _():
    import sys

    # WASM has no site-packages, so divtel has to come from a wheel.
    if sys.platform == "emscripten":
        import json

        import micropip
        from pyodide.http import pyfetch

        # Pyodide runs from a worker under assets/, so step up one level
        # to reach the site root where the wheel is published.
        base = "../pypi/"
        # The wheel's filename carries a version, so read it from a
        # manifest instead of hardcoding it.
        response = await pyfetch(base + "manifest.json")
        if response.status != 200:
            raise RuntimeError(
                f"could not fetch {base}manifest.json (HTTP {response.status}). "
                "The divtel wheel is missing from the published site."
            )
        manifest = json.loads(await response.string())
        await micropip.install(base + manifest["divtel"])
    return

@app.cell(hide_code=True)
def _():
    from importlib.resources import files

    import astropy.units as u
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.coordinates import AltAz, SkyCoord, get_body
    from astropy.utils import iers

    from divtel.layout import load_array
    from divtel.observation import Observation, pointing_coord
    from divtel.visualization import display_hyper_fov

    # Pyodide has no network, so skip astropy's default fetch of the current
    # Earth-orientation table. The bundled table is accurate to well under
    # an arcsecond near the release date, plenty for pointing an array.
    iers.conf.auto_download = False

    plt.rcParams["savefig.format"] = "svg"
    return (
        AltAz, Observation, SkyCoord, display_hyper_fov, files, get_body,
        load_array, mo, np, plt, pointing_coord, u,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.Html(
        """<style>
          img, svg { max-width: 100%; height: auto; }
        </style>"""
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Observing a real source

        The geometry so far lives in a frame fixed to the ground: an array
        pointed at 70 degrees altitude stays there whatever the hour. A
        source doesn't. It rises, crosses the sky, and sets, so the direction
        you point in changes all night.

        `Observation` bridges the two: give it a place and a time, and it
        turns a sky position into the altitude and azimuth
        `divergent_pointing` wants.

        Pick a target and a night:
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    # Hardcoded rather than resolved by name, since this runs in a browser
    # sandbox with no network. All ICRS, in degrees.
    SOURCES = {
        "Crab Nebula": (83.633, 22.015),
        "Markarian 421": (166.114, 38.209),
        "Markarian 501": (253.468, 39.760),
        "Galactic centre": (266.417, -29.008),
        "PKS 2155-304": (329.717, -30.226),
        "Centaurus A": (201.365, -43.019),
    }

    source = mo.ui.dropdown(
        list(SOURCES), value="Crab Nebula", label="source", full_width=True
    )
    site = mo.ui.dropdown(
        {"CTAO-North (La Palma)": "north", "CTAO-South (Paranal)": "south"},
        value="CTAO-North (La Palma)", label="site", full_width=True,
    )
    night = mo.ui.dropdown(
        ["2026-03-01", "2026-06-01", "2026-09-01", "2026-12-01"],
        value="2026-03-01", label="night of", full_width=True,
    )
    mo.hstack([source, site, night], widths="equal")
    return SOURCES, night, site, source


@app.cell(hide_code=True)
def _(AltAz, Observation, SOURCES, SkyCoord, get_body, night, np, site, source, u):
    TARGET = SkyCoord(
        ra=SOURCES[source.value][0] * u.deg, dec=SOURCES[source.value][1] * u.deg
    )

    def _solar_midnight():
        """
        The middle of the night, not midnight UTC.

        La Palma sits 1.2 hours west of Greenwich, Paranal 4.7, and the
        sun's lowest point drifts by up to a further quarter hour over the
        year. Centring the plot on 00:00 UTC would push the dark band off to
        one side, so one coarse ephemeris pass finds where the sun is
        actually lowest instead.
        """
        start = Observation(site=site.value, time=f"{night.value}T00:00:00")
        coarse = start.time + np.linspace(-12, 12, 49) * u.hour
        sun = get_body("sun", coarse, location=start.location).transform_to(
            AltAz(obstime=coarse, location=start.location)
        )
        return Observation(site=site.value, time=coarse[sun.alt.argmin()])

    MIDNIGHT = _solar_midnight()
    return MIDNIGHT, TARGET


@app.cell(hide_code=True)
def _(AltAz, MIDNIGHT, TARGET, get_body, np, u):
    HOURS = np.linspace(-12, 12, 97) * u.hour
    TIMES = MIDNIGHT.time + HOURS
    FRAME = AltAz(obstime=TIMES, location=MIDNIGHT.location)

    SUN = get_body("sun", TIMES, location=MIDNIGHT.location).transform_to(FRAME)
    MOON = get_body("moon", TIMES, location=MIDNIGHT.location).transform_to(FRAME)
    TRACK = TARGET.transform_to(FRAME)

    # Astronomical twilight: the sun 18 degrees below the horizon is the
    # conventional start of usable dark time for Cherenkov astronomy.
    DARK = SUN.alt < -18 * u.deg
    return DARK, HOURS, MOON, SUN, TRACK


@app.cell(hide_code=True)
def _(DARK, HOURS, MOON, SUN, TRACK, np, plt, source, u):
    def _navigation_plot():
        hours = HOURS.to_value(u.hour)
        fig, ax = plt.subplots(figsize=(10, 5))

        # Shade twilight and full dark behind everything else.
        ax.fill_between(hours, -90, 90, where=SUN.alt < 0 * u.deg,
                        color="0.75", zorder=0)
        ax.fill_between(hours, -90, 90, where=DARK, color="0.45", zorder=0)

        ax.plot(hours, SUN.alt.to_value(u.deg), color="tab:orange", label="Sun")
        ax.plot(hours, MOON.alt.to_value(u.deg), color="0.3", ls="--", label="Moon")
        scatter = ax.scatter(hours, TRACK.alt.to_value(u.deg), c=TRACK.az.to_value(u.deg),
                             cmap="twilight", s=12, zorder=3, label=source.value)
        fig.colorbar(scatter, ax=ax, label="source azimuth [deg]")

        ax.axhline(0, color="k", lw=0.8)
        ax.set_xlim(-12, 12)
        ax.set_ylim(-90, 90)
        ax.set_xticks(np.arange(-12, 13, 3))
        ax.set_xlabel("hours from the middle of the night")
        ax.set_ylabel("altitude [deg]")
        ax.set_title("what is up, and when")
        ax.legend(loc="upper left", frameon=False)
        fig.tight_layout()
        return fig

    _navigation_plot()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        The plot centres on the middle of the night, not midnight UTC at
        either site. The dark band is astronomical night, sun more than 18
        degrees down; the light band is twilight. The source track is
        coloured by azimuth, so you can see it swing across the sky as well
        as rise and fall.

        A source is worth observing when it's high *and* the sun is down,
        and those two don't always overlap: try the Galactic centre from
        CTAO-North, which barely clears the horizon, against the same
        source from CTAO-South.

        ## What the array sees while it tracks

        Now point the array at it, with `div` held fixed. The only thing
        changing is where the source has got to.
        """
    )
    return


@app.cell(hide_code=True)
def _(files, load_array):
    ARRAY = load_array(
        files("divtel") / "data" / "cta-north-lapalma-alpha-prod6.ecsv"
    )
    return (ARRAY,)


@app.cell(hide_code=True)
def _(mo):
    div = mo.ui.slider(
        0, 0.1, step=0.005, value=0.02, label="divergence", show_value=True,
        full_width=True,
    )
    hour = mo.ui.slider(
        -6, 6, step=0.5, value=0, label="hours from the middle of the night", show_value=True,
        full_width=True,
    )
    mo.vstack([div, hour])
    return div, hour


@app.cell(hide_code=True)
def _(ARRAY, MIDNIGHT, TARGET, display_hyper_fov, div, hour, mo, plt, u):
    def _snapshot():
        moment = MIDNIGHT.after(hour.value * u.hour)
        alt, az = moment.altaz_of(TARGET)

        if alt < 0 * u.deg:
            return mo.md(
                f"### The source is below the horizon at "
                f"{moment.time.isot[11:16]} UTC (altitude {alt:.1f}).\n\n"
                "Move the time slider, or pick a source that is up on this night."
            )

        ARRAY.divergent_pointing(div.value, alt, az)
        area = ARRAY.hyper_fov()[0]
        mean, _ = ARRAY.multiplicity_moments()

        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        display_hyper_fov(ARRAY, ax=ax)
        fig.tight_layout()

        return mo.hstack([
            fig,
            mo.md(
                f"**{moment.time.isot[11:16]} UTC**\n\n"
                f"- altitude **{alt:.1f}**\n"
                f"- azimuth **{az:.1f}**\n"
                f"- hyper FoV **{area:.1f}**\n"
                f"- mean multiplicity **{mean:.2f}**\n"
            ),
        ], widths=[2, 1], align="center")

    _snapshot()
    return


@app.cell(hide_code=True)
def _(ARRAY, MIDNIGHT, TARGET, div, np, plt, source, u):
    def _across_the_night():
        hours = np.arange(-6, 6.01, 0.5)
        times, areas, means = [], [], []
        for offset in hours:
            moment = MIDNIGHT.after(offset * u.hour)
            alt, az = moment.altaz_of(TARGET)
            if alt <= 0 * u.deg:
                continue
            ARRAY.divergent_pointing(div.value, alt, az)
            times.append(offset)
            areas.append(ARRAY.hyper_fov()[0].to_value(u.deg**2))
            means.append(ARRAY.multiplicity_moments()[0])

        fig, ax = plt.subplots(figsize=(9, 4.5))
        if not times:
            ax.text(0.5, 0.5, f"{source.value} does not rise on this night",
                    ha="center", va="center", transform=ax.transAxes)
            return fig

        ax.plot(times, areas, marker="o", color="tab:blue")
        ax.set_xlabel("hours from the middle of the night")
        ax.set_ylabel("hyper FoV [deg$^2$]", color="tab:blue")
        ax.tick_params(axis="y", labelcolor="tab:blue")

        twin = ax.twinx()
        twin.plot(times, means, marker="s", color="tab:red")
        twin.set_ylabel("mean multiplicity", color="tab:red")
        twin.tick_params(axis="y", labelcolor="tab:red")

        ax.set_title(f"tracking {source.value} at div = {div.value}")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        return fig

    _across_the_night()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        The covered area isn't constant while the array tracks, even though
        `div` never changes.

        The array is **flat**, on a hillside, pointing in a **tilted**
        direction. As the source drops toward the horizon, the array's
        extent along the pointing axis grows relative to its extent across
        it, so telescope offsets from the barycenter increasingly lie
        *along* the line of sight rather than across it. Divergence is
        driven by the across-axis part, so the fan closes up, the discs pile
        back onto each other, the area shrinks, and multiplicity climbs.

        It's easy to assume a fixed `div` buys a fixed field of view. It
        doesn't. Holding a hyper FoV through an observation means adjusting
        `div` as the source moves, the question the next notebook takes up.
        """
    )
    return


if __name__ == "__main__":
    app.run()

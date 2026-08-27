import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")


@app.cell(hide_code=True)
async def _():
    import sys

    # In the WASM build there is no site-packages: Python is Pyodide, and
    # divtel has to be fetched as a wheel published next to this page. The
    # released sdist on PyPI cannot be used -- Pyodide installs wheels only.
    if sys.platform == "emscripten":
        import json

        import micropip
        from pyodide.http import pyfetch

        # Pyodide runs in a worker loaded from <export root>/assets/, so
        # relative URLs resolve against assets/ rather than the page. marimo
        # always emits index.html alongside assets/, so stepping up one level
        # reaches the export root whatever subpath the site is deployed under.
        base = "../pypi/"
        # setuptools_scm bakes the version into the wheel filename, so the
        # docs build writes a manifest rather than us hardcoding it here.
        response = await pyfetch(base + "manifest.json")
        if response.status != 200:
            # Without this the JSON decode fails on a 404 page and the reader
            # is told only "see console for details".
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
    from astropy.coordinates import SkyCoord
    from astropy.utils import iers

    from divtel.layout import load_array
    from divtel.observation import Observation

    # No network in Pyodide; fall back to the table bundled with astropy.
    iers.conf.auto_download = False

    plt.rcParams["savefig.format"] = "svg"
    return Observation, SkyCoord, files, load_array, mo, np, plt, u


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
        # Choosing div

        The guide is blunt about what `div` is: an engineering knob, not a
        physical angle. It sets where a virtual point **G** sits behind the
        array, and every telescope points radially away from it. Nothing about
        the number tells you which value to run an observation at.

        This notebook works out an answer, by turning the question round. Do
        not ask what `div` to use. Ask **what you want from the array**, and
        let that pick `div`.

        Two quantities are in tension, and divergence trades one against the
        other:

        - the **hyper FoV**, the sky the array covers
        - the **multiplicity**, how many telescopes see any given part of it

        Only patches seen by two or more telescopes can be reconstructed
        stereoscopically, so the area worth quoting is the one above that cut,
        and the multiplicity must not be allowed to collapse to one.
        """
    )
    return


@app.cell(hide_code=True)
def _(files, load_array):
    ARRAY = load_array(files("divtel") / "data" / "la_palma_4LST_15MST.ecsv")
    # Sampled finely enough to interpolate a threshold crossing, coarsely
    # enough that the whole grid is a couple of seconds of geometry.
    DIVS = [0.0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.06]
    ALTITUDES = [80, 70, 60, 50, 40, 30, 20]
    AZIMUTH = 180
    return ALTITUDES, ARRAY, AZIMUTH, DIVS


@app.cell(hide_code=True)
def _(ALTITUDES, ARRAY, AZIMUTH, DIVS, np, u):
    def _scan():
        """hFoV gain over parallel pointing, and mean multiplicity, on a grid."""
        gain = np.zeros((len(ALTITUDES), len(DIVS)))
        multiplicity = np.zeros_like(gain)
        for i, alt in enumerate(ALTITUDES):
            ARRAY.divergent_pointing(0.0, alt * u.deg, AZIMUTH * u.deg)
            parallel = ARRAY.hyper_fov(m_cut=2)[0].to_value(u.deg**2)
            for j, div in enumerate(DIVS):
                ARRAY.divergent_pointing(div, alt * u.deg, AZIMUTH * u.deg)
                area, patches = ARRAY.hyper_fov(m_cut=2)
                gain[i, j] = area.to_value(u.deg**2) / parallel
                multiplicity[i, j] = ARRAY.multiplicity_moments(patches=patches)[0]
        return gain, multiplicity

    GAIN, MULTIPLICITY = _scan()
    return GAIN, MULTIPLICITY


@app.cell(hide_code=True)
def _(ALTITUDES, DIVS, GAIN, MULTIPLICITY, plt):
    def _curves():
        fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.5))
        colors = plt.get_cmap("viridis")

        for i, alt in enumerate(ALTITUDES):
            color = colors(i / max(len(ALTITUDES) - 1, 1))
            left.plot(DIVS, GAIN[i], marker="o", color=color, label=f"{alt}°")
            right.plot(DIVS, MULTIPLICITY[i], marker="o", color=color)

        left.set_xlabel("div")
        left.set_ylabel("hyper FoV gain over parallel")
        left.set_title("what divergence buys")
        left.grid(alpha=0.3)
        left.legend(title="altitude", frameon=False, ncol=2, fontsize="small")

        right.axhline(2, color="k", ls="--", lw=1)
        right.text(DIVS[-1], 2.1, "stereo limit", ha="right", fontsize="small")
        right.set_xlabel("div")
        right.set_ylabel("mean multiplicity")
        right.set_title("what it costs")
        right.set_yscale("log")
        right.grid(alpha=0.3, which="both")

        fig.tight_layout()
        return fig

    _curves()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Both panels are ordered by altitude, and the ordering is the point.

        **The same `div` buys less at low altitude.** At `div = 0.02` the array
        covers about 2.9 times its parallel field when pointing near the
        zenith, but only about 2.3 times when pointing at 20 degrees. The array
        is flat and the sky is not: as the pointing tips over, more of each
        telescope's offset from the barycenter lies *along* the line of sight
        rather than across it, and only the across-axis part drives divergence.
        The fan closes up.

        So a fixed `div` does not hold a fixed field of view while you track a
        source. It quietly gives you less as the source sets.

        ## Turning it round

        Pick a target instead, and solve for the `div` that meets it.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    target = mo.ui.slider(
        1.5, 4.0, step=0.1, value=2.5, label="wanted hyper FoV gain",
        show_value=True, full_width=True,
    )
    floor = mo.ui.slider(
        1.5, 6.0, step=0.5, value=3.0, label="minimum mean multiplicity",
        show_value=True, full_width=True,
    )
    mo.vstack([target, floor])
    return floor, target


@app.cell(hide_code=True)
def _(ALTITUDES, DIVS, GAIN, MULTIPLICITY, np):
    def required_div(target, floor):
        """
        Smallest div reaching the wanted gain without dropping multiplicity
        below the floor. NaN where no sampled div can do both.
        """
        chosen, achieved = [], []
        for i in range(len(ALTITUDES)):
            allowed = [j for j in range(len(DIVS)) if MULTIPLICITY[i, j] >= floor]
            meets = [j for j in allowed if GAIN[i, j] >= target]
            if meets:
                j = min(meets)
                # Linear interpolation between the sampled points, so the answer
                # is not quantised to the grid.
                if j > 0:
                    g0, g1 = GAIN[i, j - 1], GAIN[i, j]
                    frac = (target - g0) / (g1 - g0) if g1 > g0 else 0.0
                    chosen.append(DIVS[j - 1] + frac * (DIVS[j] - DIVS[j - 1]))
                else:
                    chosen.append(DIVS[j])
                achieved.append(target)
            else:
                # Cannot reach the target without breaking the multiplicity
                # floor: report the best that stays legal.
                if allowed:
                    j = max(allowed)
                    chosen.append(DIVS[j])
                    achieved.append(GAIN[i, j])
                else:
                    chosen.append(np.nan)
                    achieved.append(np.nan)
        return np.array(chosen), np.array(achieved)

    return (required_div,)


@app.cell(hide_code=True)
def _(ALTITUDES, required_div, floor, np, plt, target):
    def _plot_required():
        chosen, achieved = required_div(target.value, floor.value)
        short = achieved < target.value - 1e-9

        fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.5))

        left.plot(ALTITUDES, chosen, marker="o", color="tab:blue")
        left.set_xlabel("altitude [deg]")
        left.set_ylabel("div needed")
        left.set_title(f"div to reach {target.value:.1f}x the parallel field")
        left.grid(alpha=0.3)
        left.invert_xaxis()

        right.plot(ALTITUDES, achieved, marker="o", color="tab:green", label="achieved")
        right.axhline(target.value, color="k", ls="--", lw=1, label="wanted")
        if short.any():
            right.plot(np.array(ALTITUDES)[short], achieved[short], "rx",
                       markersize=11, label="multiplicity floor binds")
        right.set_xlabel("altitude [deg]")
        right.set_ylabel("hyper FoV gain")
        right.set_title("is the target reachable?")
        right.grid(alpha=0.3)
        right.invert_xaxis()
        right.legend(frameon=False, fontsize="small")

        fig.tight_layout()
        return fig

    _plot_required()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        The left panel is the practical answer: **the `div` you need climbs as
        the source drops.** Holding a constant field of view through a night is
        not a matter of setting `div` once.

        The right panel is the honest caveat. Push the wanted gain up, or the
        multiplicity floor up, and the two constraints eventually collide:
        crosses mark altitudes where no divergence reaches the target without
        thinning the array past the floor. That collision is not a limitation
        of the software. It is the trade itself — there is only so much sky
        nineteen telescopes can cover and still see any of it twice.

        ## Across the sky

        Altitude is most of the story, but not all of it: the array is not
        azimuthally symmetric, so the answer wobbles by a few percent with
        azimuth — up to about 12% at low altitude. The map below computes each
        sky position properly rather than reading it off the altitude curve.

        This is by far the most expensive cell in the notebook -- a few
        seconds of geometry natively, and rather longer through a Python
        compiled to WebAssembly -- so it waits to be asked.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    run = mo.ui.run_button(
        label="compute the sky map (up to a minute in the browser)"
    )
    # A cell's trailing expression is what marimo renders, so this is the
    # button appearing on the page rather than a statement with no effect.
    run  # noqa: B018
    return (run,)


@app.cell(hide_code=True)
def _(ARRAY, Observation, SkyCoord, floor, np, run, target, u):
    def _sky_scan():
        if not run.value:
            return None

        obs = Observation(site="north", time="2026-03-01T23:00:00")
        ra_grid = np.linspace(0, 330, 12)
        dec_grid = np.linspace(-40, 80, 7)
        divs = [0.005, 0.01, 0.02, 0.03, 0.05]

        chosen = np.full((len(dec_grid), len(ra_grid)), np.nan)
        cost = np.full_like(chosen, np.nan)
        altitude = np.full_like(chosen, np.nan)

        for row, dec in enumerate(dec_grid):
            for col, ra in enumerate(ra_grid):
                alt, az = obs.altaz_of(SkyCoord(ra=ra * u.deg, dec=dec * u.deg))
                altitude[row, col] = alt.to_value(u.deg)
                if alt < 20 * u.deg:
                    # Below about 20 degrees the atmosphere makes the pointing
                    # academic, so there is nothing to choose.
                    continue

                ARRAY.divergent_pointing(0.0, alt, az)
                parallel_area, parallel_patches = ARRAY.hyper_fov(m_cut=2)
                parallel = parallel_area.to_value(u.deg**2)

                # Walk up the divergences until the target is met, then
                # interpolate between the bracketing pair rather than reporting
                # whichever grid value happened to clear it -- the same trick
                # the altitude curve above uses, and without it the map is
                # quantised to the handful of divs actually sampled.
                last_div = 0.0
                last_gain = 1.0
                last_mean = ARRAY.multiplicity_moments(patches=parallel_patches)[0]
                for div in divs:
                    ARRAY.divergent_pointing(div, alt, az)
                    area, patches = ARRAY.hyper_fov(m_cut=2)
                    gain = area.to_value(u.deg**2) / parallel
                    mean = ARRAY.multiplicity_moments(patches=patches)[0]

                    if mean < floor.value:
                        # This div breaks the multiplicity floor; the best legal
                        # answer is the previous one.
                        break

                    if gain >= target.value:
                        span = gain - last_gain
                        frac = (target.value - last_gain) / span if span > 0 else 0.0
                        chosen[row, col] = last_div + frac * (div - last_div)
                        cost[row, col] = last_mean + frac * (mean - last_mean)
                        break

                    last_div, last_gain, last_mean = div, gain, mean
                    chosen[row, col] = div
                    cost[row, col] = mean

        return ra_grid, dec_grid, chosen, cost, altitude

    SKY = _sky_scan()
    return (SKY,)


@app.cell(hide_code=True)
def _(SKY, mo, plt, target):
    def _sky_map():
        if SKY is None:
            return mo.md(
                "*Press the button above to run the scan.*"
            )
        ra_grid, dec_grid, chosen, cost, altitude = SKY
        extent = [ra_grid.min(), ra_grid.max(), dec_grid.min(), dec_grid.max()]

        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        for ax, data, label, cmap in (
            (axes[0], altitude, "altitude [deg]", "cividis"),
            (axes[1], chosen, "div chosen", "viridis"),
            (axes[2], cost, "mean multiplicity", "magma"),
        ):
            image = ax.imshow(data, origin="lower", extent=extent, aspect="auto", cmap=cmap)
            fig.colorbar(image, ax=ax, label=label)
            ax.set_xlabel("RA [deg]")
            ax.set_ylabel("Dec [deg]")

        axes[0].set_title("where the sky is, at 23:00 UTC")
        axes[1].set_title(f"div for a {target.value:.1f}x field")
        axes[2].set_title("multiplicity you end up with")
        fig.suptitle("CTAO-North, 2026-03-01 23:00 UTC — blank where the source is below 20°")
        fig.tight_layout()
        return fig

    _sky_map()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
The first two panels line up: wherever the sky is low, the map asks
        for more divergence, exactly as the altitude curve said it would.

        The third panel is the payoff, and it is the reason to bother adapting
        `div` at all. Across the whole visible sky the multiplicity you end up
        with barely moves. Holding the field of view fixed by raising `div` as
        the source drops turns out to hold the *multiplicity* fixed too: over
        altitudes from 80 down to 20 degrees, `div` has to climb about 40%,
        from 0.016 to 0.022, while the mean multiplicity drifts only from 5.3
        to 5.6.

        That is worth stating as the practical result. **An adaptive `div` gives
        you a stable array all night — same coverage, same depth — where a
        fixed one gives you neither.**

        The blank band is sky that has not risen. At a fixed moment more than
        half of this map is simply unavailable, which is the other half of
        scheduling an observation.

        ## What to take away

        1. **`div` is not a setting, it is the answer to a question.** Decide
           what field of view you want and what multiplicity you can live with,
           and those two fix it.
        2. **The answer moves while you observe.** A source at 70 degrees and
           the same source at 25 degrees need different divergence for the same
           coverage. Holding `div` fixed through a night silently gives away
           field of view as the source sets.
        3. **Adapting it costs nothing in depth.** Raising `div` to hold the
           field of view as a source sets keeps the multiplicity almost
           constant too. You are not trading depth for stability; you are
           getting both, which a fixed `div` does not give you.
        4. **The two constraints can still be incompatible.** When the crosses
           appear, no divergence satisfies both; the array is not big enough to
           cover that much sky stereoscopically, and something has to give.

        The numbers here are specific to the 4 LST + 15 MST La Palma layout and
        to a stereo cut of two. Change the layout, the cut or the floor and the
        curves move — but the shape of the argument does not.
        """
    )
    return


if __name__ == "__main__":
    app.run()

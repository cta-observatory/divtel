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

    from divtel import strategy, visualization
    from divtel.layout import load_array
    from divtel.pointing import best_pointing
    from divtel.region import SkyRegion

    # SVG scales to its container; marimo's PNG path stamps a fixed pixel
    # width that overflows a frame narrower than the figure.
    plt.rcParams["savefig.format"] = "svg"
    plt.rcParams["font.size"] = 9

    DATA = files("divtel") / "data"
    return (
        DATA,
        SkyRegion,
        best_pointing,
        load_array,
        mo,
        np,
        plt,
        strategy,
        u,
        visualization,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.Html("""<style> img, svg { max-width: 100%; height: auto; } </style>""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # Covering GW170817

    A gravitational-wave localization is a probability distribution, not a
    point, and the early ones are **arcs**: long, thin, and lumpy. There are
    four ways to put an array on one, and they are not equally good.

    Pick a map, a site and a strategy, and watch the array's *multiplicity* —
    how many cameras have each direction in frame — change shape. Two or more
    telescopes can reconstruct a shower; one cannot, so grey is blind and the
    palest blue is nearly so.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    MAPS = {
        "BAYESTAR H+L — the first alert, +0.7 h": "bayestar_hl",
        "BAYESTAR H+L+V — with Virgo, +5.2 h": "bayestar_hlv",
        "LALInference (preliminary), +16.3 h": "lalinference_prelim",
    }

    which_map = mo.ui.dropdown(MAPS, value=list(MAPS)[0], label="localization")
    level = mo.ui.dropdown({"50 %": 50, "90 %": 90, "95 %": 95}, value="90 %",
                           label="credible level")
    site = mo.ui.dropdown({"CTAO-South (14 MST + 37 SST)": "south",
                           "CTAO-North (4 LST + 9 MST)": "north"},
                          value="CTAO-South (14 MST + 37 SST)", label="array")
    altitude = mo.ui.slider(20, 85, value=60, step=5, label="placed at altitude",
                            show_value=True)

    mo.hstack([which_map, level, site, altitude], justify="start", gap=1.5)
    return altitude, level, site, which_map


@app.cell(hide_code=True)
def _(DATA, SkyRegion, altitude, level, load_array, site, u, which_map):
    LAYOUTS = {"north": "cta-north-lapalma-alpha-prod6.ecsv",
               "south": "cta-south-paranal-alpha-prod6.ecsv"}

    array = load_array(DATA / LAYOUTS[site.value])

    # The region is rotated rigidly to the chosen altitude, then cut at the
    # horizon: an arc a hundred degrees long can run underground, and an
    # optimiser handed that will aim telescopes there and report the coverage.
    _whole = SkyRegion.from_table(
        DATA / "gw170817" / f"{which_map.value}_{level.value}.ecsv.gz")
    region, underground = _whole.place(altitude.value * u.deg, 180 * u.deg).visible_part()
    return array, region, underground


@app.cell(hide_code=True)
def _(mo):
    STRATEGIES = {
        "parallel — everything on one spot": "parallel",
        "divergent — one number, fanned out": "divergent",
        "tiling — one array, moved in steps": "tiling",
        "sub-arrays — split it, point the pieces": "subarrays",
        "sub-arrays, sized by probability": "weighted",
        "shaped — every telescope on its own": "shaped",
    }

    how = mo.ui.dropdown(STRATEGIES, value=list(STRATEGIES)[3],
                         label="pointing strategy")
    div = mo.ui.slider(0.0, 0.30, value=0.03, step=0.005, label="divergence d",
                       show_value=True)
    groups = mo.ui.slider(2, 17, value=8, step=1,
                          label="tiles, or sub-arrays", show_value=True)
    beta = mo.ui.slider(0.0, 2.0, value=1.0, step=0.25,
                        label="β — how hard to track the probability",
                        show_value=True)

    mo.hstack([how, div, groups, beta], justify="start", gap=1.5)
    return beta, div, groups, how


@app.cell(hide_code=True)
def _(mo, how):
    # Shaped pointing is a coordinate ascent over every telescope at once. In a
    # browser that is seconds rather than milliseconds, so it does not run on a
    # slider drag -- and the settings below are cut down from the ones the
    # study page's figures were computed with.
    solve = mo.ui.run_button(label="solve the shaped pointing")
    solve if how.value == "shaped" else mo.md("")
    return (solve,)


@app.cell
def _(array, best_pointing, beta, div, groups, how, region, solve, strategy):
    exposures, note = 1, ""

    if how.value == "parallel":
        aim = best_pointing(array, region, 0.0)
        note = f"aimed at alt {aim['alt']:.1f}, az {aim['az']:.1f}"

    elif how.value == "divergent":
        best_pointing(array, region, div.value)

    elif how.value == "tiling":
        tiles = strategy.tile_region(array, region, pointings=groups.value)
        # The picture shows one tile; the coverage below is what the whole
        # sequence reaches, which costs one exposure per tile.
        array.divergent_pointing(0.0, tiles[0]["alt"], tiles[0]["az"])
        exposures = len(tiles)
        note = (f"{100 * tiles[-1]['cumulative']:.0f} % after all "
                f"{len(tiles)} tiles; the panel shows the first")

    elif how.value == "subarrays":
        summary = strategy.point_subarrays(array, region, groups.value)
        note = f"group sizes {summary['sizes']}"

    elif how.value == "weighted":
        summary = strategy.weighted_split(array, region, groups.value)
        note = f"group sizes {summary['sizes']}"

    elif how.value == "shaped":
        if solve.value:
            # Cut down for the browser: fewer candidate pointings, fewer scored
            # directions, fewer sweeps and one warm start instead of seven.
            strategy.shaped_pointing(array, region, beta=beta.value,
                                     candidates=400, search_pixels=600,
                                     sweeps=4, rounds=1, warm_starts=(6,))
            note = "cut-down solve; the study page uses the full budget"
        else:
            best_pointing(array, region, 0.0)
            note = "press the button to solve"

    scored = strategy.describe(array, region, m_cut=2)
    profile = strategy.multiplicity_by_probability(array, region)

    # Every strategy above re-points `array` in place, and marimo re-runs a cell
    # when a variable it reads is *redefined* -- mutation is invisible to it. So
    # hand the pointed array on under a new name, or the figures below go stale
    # while the numbers beside them update.
    pointed = array
    return exposures, note, pointed, profile, scored


@app.cell(hide_code=True)
def _(exposures, mo, np, pointed, region, scored, underground, u):
    # How many telescopes have any part of the region in frame. It is the
    # number that separates the strategies most sharply: divergence spreads
    # the array so wide that a third of it ends up staring at empty sky.
    _cos_radii = np.cos([t.fov_radius.to_value(u.rad)
                         for t in pointed.telescopes])
    _on_target = int(((region.directions @ pointed.pointing_vectors.T)
                      >= _cos_radii).any(axis=0).sum())

    mo.hstack([
        mo.stat(label="covered in stereo",
                value=f"{100 * scored['covered_stereo']:.0f} %",
                caption=f"of the region above the horizon "
                        f"({100 * underground:.0f} % of the map is not)"),
        mo.stat(label="telescopes per shower",
                value=f"{scored['mean_covered']:.1f}",
                caption="averaged over the part covered"),
        mo.stat(label="telescopes on target",
                value=f"{_on_target} of {len(pointed.telescopes)}",
                caption="the rest are pointed at empty sky"),
        mo.stat(label="exposures", value=str(exposures),
                caption="the window is divided between them"),
    ], justify="start", gap=2)
    return


@app.cell(hide_code=True)
def _(mo, note):
    mo.md(f"*{note}*") if note else mo.md("")
    return


@app.cell(hide_code=True)
def _(plt, pointed, region, visualization):
    # Where the cameras actually are, over the probability they are covering.
    # The panel below shades the result; this one shows the cause.
    _frame = visualization.projection_frame(region)
    _figure, _ax = plt.subplots(figsize=(11.0, 4.0))

    visualization.probability_over_region(region, ax=_ax, frame=_frame, size=3.0)
    _rims = visualization.camera_rims(pointed, ax=_ax, frame=_frame,
                                      colors="type",
                                      type_names=("SST", "MST"))
    visualization.frame_on(
        _ax, _rims, visualization.region_extent(region, frame=_frame))

    _legend = _ax.legend(frameon=False, fontsize=8, loc="upper left", ncol=2)
    for _handle in _legend.legend_handles:
        if hasattr(_handle, "set_sizes"):
            _handle.set_sizes([26])
    _ax.set_title("the cameras, over the probability they cover", loc="left")
    _figure.tight_layout()
    _figure
    return


@app.cell(hide_code=True)
def _(plt, pointed, region, visualization):
    # Its own figure, so the equal-aspect panel can take the full width rather
    # than sharing it with a bar chart that wants a different shape.
    _figure, _ax = plt.subplots(
        figsize=(11.0, 1.6 + 11.0 * visualization.region_span(region)))
    visualization.multiplicity_over_region(pointed, region, ax=_ax, size=6.0)
    _ax.set_title("what that delivers: telescopes per direction", loc="left")
    _figure.tight_layout()
    _figure
    return


@app.cell(hide_code=True)
def _(plt, profile, visualization):
    _figure, _ax = plt.subplots(figsize=(11.0, 3.0))
    visualization.multiplicity_by_probability({"this configuration": profile},
                                              ax=_ax)
    _figure.tight_layout()
    _figure
    return


@app.cell(hide_code=True)
def _(mo, profile, scored):
    mo.md(
        f"""
    **Reading the three panels.** The first draws every telescope's field of
    view over the localization, shaded darkest where the source most likely is.
    Rims are coloured by camera size, so you can watch free pointing send the
    narrow cameras where the map is sharp — something neither divergence nor an
    even split can do. All three panels are equal-area, so a wide patch of
    shallow coverage looks as big as it is.

    The second shades the map by how many cameras land on each direction; grey
    is sky no camera reaches. The third cuts the map into fifths of
    *probability* — the last fifth is the small bright core, the first is the
    wide faint skirt — and asks how deeply each is watched.

    A flat profile is an array spread evenly over the region, which is only the
    right answer if the map is flat. This one runs
    **{profile[0]:.1f} → {profile[-1]:.1f}**, and its mismatch against the
    probability is **{scored['mismatch']:.3f}** (zero would be depth exactly
    proportional to probability over the part that clears the stereo cut).

    Read that number beside the coverage and never alone: an array that
    abandons everything but the core matches the shape of what it kept,
    perfectly, and is useless.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ---

    Everything here is geometry: where cameras point and what they cover. There
    is no effective area, no energy threshold and no background, so none of
    these percentages is a sensitivity. Coverage is necessary for a detection
    and nowhere near sufficient.

    The full argument, with the numbers, is in
    [the study page](../../gw170817.html).
    """
    )
    return


if __name__ == "__main__":
    app.run()

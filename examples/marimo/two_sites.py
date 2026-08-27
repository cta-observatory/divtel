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

    from divtel.layout import load_array
    from divtel.visualization import display_hyper_fov

    # SVG scales to its container; marimo's PNG path stamps a fixed pixel
    # width that overflows a frame narrower than the figure.
    plt.rcParams["savefig.format"] = "svg"
    return display_hyper_fov, files, load_array, mo, plt, u


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
    mo.md(r"""
    # Two sites

    CTAO has two sites, one in each hemisphere. The North site, on La Palma, is four
    Large-Sized Telescopes surrounded by nine Medium-Sized ones. The
    South site, at Paranal, has fourteen MSTs and
    thirty-seven Small-Sized Telescopes, the SSTs picking up the
    high-energy end of the spectrum that the North doesn't chase.
    """)
    return


@app.cell(hide_code=True)
def _(files, load_array):
    NORTH = load_array(
        files("divtel") / "data" / "cta-north-lapalma-alpha-prod6.ecsv"
    )
    SOUTH = load_array(
        files("divtel") / "data" / "cta-south-paranal-alpha-prod6.ecsv"
    )
    # Telescope ids restart at 1 on each site. It's a per-file convention,
    # not one id space shared across CTAO. Group across sites by camera
    # radius, not by id, if the two ever need to be merged.
    return NORTH, SOUTH


@app.cell(hide_code=True)
def _(NORTH, SOUTH, mo, u):
    def _composition():
        rows = []
        for site_name, array in [("North (La Palma)", NORTH), ("South (Paranal)", SOUTH)]:
            radii = sorted({round(tel.fov_radius.to_value(u.deg), 2) for tel in array.telescopes})
            rows.append(
                f"| {site_name} | {len(array.telescopes)} | "
                f"{', '.join(f'{r}°' for r in radii)} |"
            )
        return mo.md(
            "| site | telescopes | camera FoV radii |\n"
            "|---|---|---|\n" + "\n".join(rows)
        )

    _composition()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The South array is both bigger and wider-eyed: its SSTs subtend
    nearly 11 degrees each, almost three times an MST camera and five
    times an LST's, since the SSTs are built to catch the biggest,
    highest-energy showers over the widest possible area rather than the
    faint, low-energy ones an LST resolves.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    div = mo.ui.slider(
        0, 0.1, step=0.005, value=0.02, label="divergence", show_value=True,
        full_width=True,
    )
    alt = mo.ui.slider(
        20, 90, step=1, value=70, label="altitude [deg]", show_value=True,
        full_width=True,
    )
    mo.vstack([div, alt])
    return alt, div


@app.cell(hide_code=True)
def _(NORTH, SOUTH, alt, div, u):
    NORTH.divergent_pointing(div.value, alt.value * u.deg, 180 * u.deg)
    SOUTH.divergent_pointing(div.value, alt.value * u.deg, 180 * u.deg)
    # A new binding, not just the mutated arrays, so cells reading POINTED
    # rerun on every slider move: marimo reacts to which names a cell
    # defines, and re-pointing NORTH/SOUTH in place doesn't define anything.
    POINTED = (div.value, alt.value)
    return (POINTED,)


@app.cell(hide_code=True)
def _(NORTH, POINTED, SOUTH, display_hyper_fov, plt):
    # POINTED isn't used directly; it's a dependency so this cell reruns
    # whenever the sliders move NORTH's and SOUTH's pointing.
    _ = POINTED

    def _both_arrays():
        fig, axes = plt.subplots(2, 2, figsize=(11, 10))
        (north_ground, south_ground), (north_sky, south_sky) = axes

        NORTH.display_2d(projection="xy", ax=north_ground)
        SOUTH.display_2d(projection="xy", ax=south_ground)
        north_ground.set_title("North, on the ground")
        south_ground.set_title("South, on the ground")

        display_hyper_fov(NORTH, ax=north_sky)
        display_hyper_fov(SOUTH, ax=south_sky)
        north_sky.set_title("North: " + north_sky.get_title())
        south_sky.set_title("South: " + south_sky.get_title())

        fig.tight_layout()
        return fig

    _both_arrays()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    South's hyper FoV runs roughly ten times North's at the same `div`
    and altitude. Four times the telescope count, cameras up to five
    times wider, and the two effects compound rather than add: an order
    of magnitude more sky.

    South also starts from a higher mean
    multiplicity at the same `div`. Stack discs from 51 telescopes
    instead of 13 and any given patch of covered sky tends to be seen by
    more of them, so South can push `div` further before multiplicity
    thins to the stereo limit. Matching the two arrays at one shared
    `div` value, the way the sliders above do, is a start, not the last
    word. The "Choosing div" tutorial picks `div` from a wanted field of
    view and a multiplicity floor instead, and that's closer to how
    North and South would each get tuned in practice.
    """)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()

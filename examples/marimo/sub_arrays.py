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
    from divtel.visualization import display_groups, display_hyper_fov

    # SVG scales to its container; marimo's PNG path stamps a fixed pixel
    # width that overflows a frame narrower than the figure.
    plt.rcParams["savefig.format"] = "svg"
    return display_groups, display_hyper_fov, files, load_array, mo, plt, u


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
        # One site, two instruments

        The CTAO North layout at La Palma is really two arrays sharing a
        field: four Large-Sized Telescopes sitting inside nine Medium-Sized
        ones, with different optics and different fields of view. An LST
        camera subtends about 2.15 degrees on the sky, an MST camera about
        3.84.

        Divergence is a *geometric* construction: every telescope points
        away from a virtual point behind the array, so how far it swings
        depends on how far it sits from the centre. The LSTs cluster near
        the middle; the MSTs spread out around them.

        This notebook is about what that asymmetry does.
        """
    )
    return


@app.cell(hide_code=True)
def _(files, load_array):
    ARRAY = load_array(
        files("divtel") / "data" / "cta-north-lapalma-alpha-prod6.ecsv"
    )
    # Ids follow the CTAO convention, so type maps directly to id range.
    # `ARRAY.group_by("camera_radius")` finds the same split unprompted,
    # since each type shares a camera.
    TYPES = {"LST": range(1, 5), "MST": range(5, 14)}
    return ARRAY, TYPES


@app.cell(hide_code=True)
def _(mo):
    div = mo.ui.slider(
        0, 0.3, step=0.005, value=0.02, label="divergence", show_value=True,
        full_width=True,
    )
    alt = mo.ui.slider(
        20, 90, step=1, value=70, label="altitude [deg]", show_value=True,
        full_width=True,
    )
    mo.vstack([div, alt])
    return alt, div


@app.cell(hide_code=True)
def _(ARRAY, TYPES, alt, div, u):
    ARRAY.divergent_pointing(div.value, alt.value * u.deg, 180 * u.deg)
    GROUPS = ARRAY.group_by(TYPES)
    return (GROUPS,)


@app.cell(hide_code=True)
def _(ARRAY, GROUPS, display_groups, display_hyper_fov, plt, u):
    def _both_views():
        fig, (ground, sky) = plt.subplots(1, 2, figsize=(11, 5))
        display_groups(GROUPS, ax=ground)
        display_hyper_fov(ARRAY, ax=sky)
        ground.set_title("on the ground")
        fig.tight_layout()
        return fig

    _both_views()
    return


@app.cell(hide_code=True)
def _(ARRAY, GROUPS, mo, u):
    def _summary():
        rows = []
        for name, group in list(GROUPS.items()) + [("both", ARRAY)]:
            area = group.hyper_fov()[0].to_value(u.deg**2)
            mean, _ = group.multiplicity_moments()
            rows.append(
                f"| {name} | {len(group.telescopes)} | {area:.1f} | {mean:.2f} |"
            )
        return mo.md(
            "| | telescopes | hyper FoV [deg²] | mean multiplicity |\n"
            "|---|---|---|---|\n" + "\n".join(rows)
        )

    _summary()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Two things to notice

        **At low divergence, the MST row and the "both" row report the same
        area.** Not a coincidence: an MST camera is wider than an LST one,
        and the two types sit on top of each other, so the sky the LSTs see
        is *inside* the sky the MSTs see. Adding the LSTs back buys no new
        sky, just multiplicity, the "both" row's mean is higher because four
        more telescopes are piled onto a patch the MSTs already watch.

        Worth knowing before quoting a hyper FoV for a mixed array: the
        number is set by the widest camera. Narrow-camera telescopes add
        depth, not width, and area alone won't show that.

        **One `div` doesn't mean one angle.** A telescope's divergence angle
        is set by how far it sits from the array centre, across the pointing
        axis:

        $$\alpha_i = \arctan\frac{|r_{\perp,i}|}{\text{norm} + r_{\parallel,i}}$$

        The LSTs cluster near the centre; the MSTs run out past 300 m. At
        `div = 0.02` the LSTs swing between 0.5 and 1.4 degrees off the mean
        pointing, the MSTs between 0.5 and 3.7, roughly two and a half times
        the angle on the same knob.

        The consequence is easy to guess wrong. It does *not* follow that
        the MSTs lose stereo overlap first: their camera is also the wider
        one, 3.84 degrees of radius against 2.15. The bigger swing and the
        bigger camera partly cancel out, so both types thin to multiplicity
        one at a similar `div`, around 0.08 to 0.1 here, even though the
        angular spread stays roughly two and a half times apart the whole
        way.
        """
    )
    return


@app.cell(hide_code=True)
def _(ARRAY, TYPES, plt, u):
    def _spread():
        divs = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3]
        curves = {name: [] for name in TYPES}
        for value in divs:
            ARRAY.divergent_pointing(value, 70 * u.deg, 180 * u.deg)
            for name, group in ARRAY.group_by(TYPES).items():
                curves[name].append(group.hyper_fov()[0].to_value(u.deg**2))

        fig, ax = plt.subplots(figsize=(7, 4.5))
        for name, areas in curves.items():
            ax.plot(divs, areas, marker="o", label=name)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("div")
        ax.set_ylabel("hyper FoV [deg$^2$]")
        ax.set_title("the two types respond to div differently")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(frameon=False)
        fig.tight_layout()
        return fig

    _spread()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        The two curves never meet. As `div` goes to zero, each type collapses
        onto a single camera, so the gap bottoms out at the ratio of the two
        camera areas, about 3.2. The MST curve then pulls away, reaching
        roughly seven times the LST area by `div = 0.3`, and both flatten
        once every telescope sees its own patch of sky, no more area to win.

        To make the two types diverge by comparable *angles* instead of a
        shared `div`, group the array and point each group on its own. Each
        sub-array is a full `Array`, with its own barycenter and its own
        `divergent_pointing`.
        """
    )
    return


if __name__ == "__main__":
    app.run()

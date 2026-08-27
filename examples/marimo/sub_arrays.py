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

    from divtel.layout import load_array
    from divtel.visualization import display_groups, display_hyper_fov

    # Render figures as SVG. marimo's PNG path stamps an explicit pixel width
    # on the image -- figsize x 100 -- which overflows a frame narrower than
    # the figure. SVG carries no such width, so it scales to its container.
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

        The CTA La Palma layout is not one array but two sharing a field: four
        Large-Sized Telescopes sitting inside fifteen Medium-Sized ones. They
        have different optics, and so different fields of view -- an LST camera
        subtends about 2.15 degrees on the sky, an MST camera about 3.84.

        Divergence acts on them differently, because it is a *geometric*
        construction: every telescope points away from a virtual point behind
        the array, so how far a telescope swings depends on how far it sits
        from the centre. The LSTs are clustered in the middle. The MSTs are
        spread out around them.

        This notebook is about what that asymmetry does.
        """
    )
    return


@app.cell(hide_code=True)
def _(files, load_array):
    ARRAY = load_array(files("divtel") / "data" / "la_palma_4LST_15MST.ecsv")
    # Ids follow the CTAO convention, so the telescope types are exactly the
    # id ranges. `ARRAY.group_by("camera_radius")` finds the same split without
    # being told it, since telescopes of a type share a camera.
    TYPES = {"LST": range(1, 5), "MST": range(5, 20)}
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
        area.** That is not a rounding coincidence: an MST camera is wider than
        an LST one and the two types sit on top of each other, so the sky the
        LSTs see is *inside* the sky the MSTs see. Adding the LSTs back buys no
        new sky at all. What it buys is multiplicity -- the "both" row's mean
        is higher than either type alone, because the LSTs pile four more
        telescopes onto a patch the MSTs were already watching.

        This is worth knowing before quoting a hyper FoV for a mixed array. The
        number is set by the widest camera; the narrow-camera telescopes are
        contributing depth, not width, and an area figure alone will not show
        that.

        **One `div` does not mean one angle.** A telescope's divergence angle
        is set by how far it sits from the array centre, across the pointing
        axis:

        $$\alpha_i = \arctan\frac{|r_{\perp,i}|}{\text{norm} + r_{\parallel,i}}$$

        The LSTs are clustered near the centre; the MSTs run out to nearly
        400 m. At `div = 0.02` the LSTs swing between 0.4 and 1.8 degrees off
        the mean pointing while the MSTs swing between 0.8 and 4.8 -- the same
        knob, roughly three times the angle.

        What is easy to guess wrong is the consequence. It does *not* follow
        that the MSTs lose their stereoscopic overlap first, because the MST
        camera is also the wider one: 3.84 degrees of radius against 2.15. The
        bigger swing and the bigger camera pull against each other and largely
        cancel, so both types thin out to multiplicity one at a broadly similar
        `div` -- around 0.05 to 0.1 on this layout. The angular spreads stay a
        factor of three apart the whole way.
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
        The two curves never meet. As `div` goes to zero each type collapses
        onto a single camera, so the gap bottoms out at the ratio of the two
        camera areas, about 3.2. From there the MST curve pulls away, reaching
        roughly twelve times the LST area by `div = 0.1`, and then both
        flatten: once every telescope of a type sees its own patch of sky,
        there is no more area to win.

        If you want the two types to diverge by comparable *angles* rather than
        under a comparable `div`, group the array and point each group on its
        own. Each sub-array is a full `Array`, with its own barycenter and its
        own `divergent_pointing`.
        """
    )
    return


if __name__ == "__main__":
    app.run()

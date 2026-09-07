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
    mo.md(r"""
    # Two sites, two instruments each

    On each site, CTAO is composed of two subarrays of different
    telescope types sharing a field: a smaller, inner group and a
    larger, outer one with a wider camera. La Palma pairs LSTs with
    MSTs; Paranal pairs MSTs with SSTs.

    Divergence is a *geometric* construction: every telescope points
    away from a virtual point behind the array, so how far it swings
    depends on how far it sits from the centre. The inner group
    clusters near the middle; the outer one spreads out around it.

    This notebook looks at what that asymmetry does, site by site.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## North: La Palma, LST + MST

    Four Large-Sized Telescopes sit inside nine Medium-Sized ones, with
    different optics and different fields of view. An LST camera
    subtends about 2.15 degrees on the sky, an MST camera about 3.84,
    so the outer, more numerous subarray also carries the wider camera.
    """)
    return


@app.cell(hide_code=True)
def _(files, load_array):
    ARRAY_N = load_array(
        files("divtel") / "data" / "cta-north-lapalma-alpha-prod6.ecsv"
    )
    # Ids follow the CTAO convention, so type maps directly to id range.
    # `ARRAY_N.group_by("camera_radius")` finds the same split unprompted,
    # since each type shares a camera.
    TYPES_N = {"LST": range(1, 5), "MST": range(5, 14)}
    return ARRAY_N, TYPES_N


@app.cell(hide_code=True)
def _(mo):
    div_n = mo.ui.slider(
        0, 0.3, step=0.005, value=0.02, label="divergence", show_value=True,
        full_width=True,
    )
    alt_n = mo.ui.slider(
        20, 90, step=1, value=70, label="altitude [deg]", show_value=True,
        full_width=True,
    )
    mo.vstack([div_n, alt_n])
    return alt_n, div_n


@app.cell(hide_code=True)
def _(ARRAY_N, TYPES_N, alt_n, div_n, u):
    ARRAY_N.divergent_pointing(div_n.value, alt_n.value * u.deg, 180 * u.deg)
    GROUPS_N = ARRAY_N.group_by(TYPES_N)
    return (GROUPS_N,)


@app.cell(hide_code=True)
def _(ARRAY_N, GROUPS_N, display_groups, display_hyper_fov, plt):
    def _both_views():
        fig, (ground, sky) = plt.subplots(1, 2, figsize=(11, 5))
        display_groups(GROUPS_N, ax=ground)
        display_hyper_fov(ARRAY_N, ax=sky)
        ground.set_title("on the ground")
        fig.tight_layout()
        return fig

    _both_views()
    return


@app.cell(hide_code=True)
def _(ARRAY_N, GROUPS_N, mo, u):
    def _summary():
        rows = []
        for name, group in list(GROUPS_N.items()) + [("both", ARRAY_N)]:
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
    mo.md(r"""
    ### Two things to notice

    **Hyper FoV here means stereo FoV.** `hyper_fov` only counts sky seen
    by two or more telescopes by default, since a shower only one
    telescope sees can't be reconstructed stereoscopically. That is why
    every row in the table above rises with `div` at first, then falls
    back towards zero: past some divergence, overlap runs out.

    **Mixing types buys back overlap the wider camera alone can't.** At
    low `div` the "both" row simply tracks the MSTs, the wider camera,
    same as before. But push `div` past the point where the MSTs stop
    overlapping *each other* and the picture changes: an MST and an LST
    can still share a patch of sky even when no two MSTs do. At
    `div = 0.1` the MSTs alone are down to 4.6 deg², yet the full array
    still covers 19.5 deg², over four times as much, entirely from
    LST-MST pairs. At `div = 0.12` the MSTs alone are essentially gone
    (0.7 deg²) while the full array still holds 5.1 deg². Cross-type
    stereo is what keeps the array useful in that gap.

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
    """)
    return


@app.cell(hide_code=True)
def _(ARRAY_N, TYPES_N, plt, u):
    def _spread():
        divs = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.12, 0.15]
        curves = {name: [] for name in list(TYPES_N) + ["both"]}
        for value in divs:
            ARRAY_N.divergent_pointing(value, 70 * u.deg, 180 * u.deg)
            for name, group in ARRAY_N.group_by(TYPES_N).items():
                curves[name].append(group.hyper_fov()[0].to_value(u.deg**2))
            curves["both"].append(ARRAY_N.hyper_fov()[0].to_value(u.deg**2))

        # Linear, not log: the stereo area actually reaches zero, which a
        # log axis can't show.
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for name, areas in curves.items():
            ax.plot(divs, areas, marker="o", label=name)
        ax.set_xlabel("div")
        ax.set_ylabel("stereo hyper FoV [deg$^2$]")
        ax.set_title("overlap rises, then runs out")
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=False)
        fig.tight_layout()
        return fig

    _spread()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Every curve rises, peaks, then falls back to zero, and they peak in
    the order you'd expect from how tightly each type is packed: the
    LSTs, closest together, peak first and lowest, around 17 deg² near
    `div = 0.02`; the MSTs peak later and higher, around 127 deg² near
    `div = 0.04`. The "both" curve tracks the MSTs closely up to their
    peak, then pulls ahead of them, held up by cross-type pairs alone
    while the MST-only curve keeps falling. By `div = 0.15` every curve
    is at zero: no two telescopes anywhere in the array still share a
    patch of sky.

    To make the two types diverge by comparable *angles* instead of a
    shared `div`, group the array and point each group on its own. Each
    sub-array is a full `Array`, with its own barycenter and its own
    `divergent_pointing`.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## South: Paranal, MST + SST

    Paranal pairs the same idea with a different pair of instruments:
    14 Medium-Sized Telescopes sit inside 37 Small-Sized ones. The SST
    mirror is by far the smaller of the two, but its focal length is
    also much shorter, so its camera ends up covering the wider patch
    of sky: about 4.4 degrees of radius against 3.75 for the MST. Same
    shape as La Palma — the outer, more numerous subarray carries the
    wider camera — but the two cameras here are much closer in size.
    """)
    return


@app.cell(hide_code=True)
def _(files, load_array):
    ARRAY_S = load_array(
        files("divtel") / "data" / "cta-south-paranal-alpha-prod6.ecsv"
    )
    TYPES_S = {"MST": range(1, 15), "SST": range(15, 52)}
    return ARRAY_S, TYPES_S


@app.cell(hide_code=True)
def _(mo):
    div_s = mo.ui.slider(
        0, 0.3, step=0.005, value=0.02, label="divergence", show_value=True,
        full_width=True,
    )
    alt_s = mo.ui.slider(
        20, 90, step=1, value=70, label="altitude [deg]", show_value=True,
        full_width=True,
    )
    mo.vstack([div_s, alt_s])
    return alt_s, div_s


@app.cell(hide_code=True)
def _(ARRAY_S, TYPES_S, alt_s, div_s, u):
    ARRAY_S.divergent_pointing(div_s.value, alt_s.value * u.deg, 180 * u.deg)
    GROUPS_S = ARRAY_S.group_by(TYPES_S)
    return (GROUPS_S,)


@app.cell(hide_code=True)
def _(ARRAY_S, GROUPS_S, display_groups, display_hyper_fov, plt):
    def _both_views():
        fig, (ground, sky) = plt.subplots(1, 2, figsize=(11, 5))
        display_groups(GROUPS_S, ax=ground)
        display_hyper_fov(ARRAY_S, ax=sky)
        ground.set_title("on the ground")
        fig.tight_layout()
        return fig

    _both_views()
    return


@app.cell(hide_code=True)
def _(ARRAY_S, GROUPS_S, mo, u):
    def _summary():
        rows = []
        for name, group in list(GROUPS_S.items()) + [("both", ARRAY_S)]:
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
    mo.md(r"""
    ### Cross-type pairs matter even more here

    At `div = 0.02` the "both" row still matches the SST row almost
    exactly, same reason as at La Palma: the MSTs' narrower sky sits
    inside the SSTs' wider one, so folding the MSTs in adds
    multiplicity, not area, while both curves are still climbing.

    The angles tell a similar story. At `div = 0.02` the MSTs swing
    between 0.2 and 3.7 degrees off the mean pointing, the SSTs between
    1.9 and 11.1 — roughly three times the spread, echoing the LST/MST
    split at La Palma even though the telescopes are different.

    Push `div` further and cross-type pairs carry far more of the
    coverage here than at La Palma. At `div = 0.1` the MSTs no longer
    overlap each other at all, and the SSTs are down to 1.7 deg², yet
    the full array still covers 65 deg² — essentially all of it from
    MST-SST pairs. At `div = 0.12` neither type overlaps itself at all,
    and the mixed array still holds 26 deg²: the entire remaining
    stereo field of view comes from telescopes of different types
    seeing the same patch. With 51 telescopes of two very differently
    sized cameras to pair up, Paranal keeps a usable stereo field far
    past the point where either subarray alone has given out.
    """)
    return


if __name__ == "__main__":
    app.run()

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


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
        base = "../wheels/"
        # setuptools_scm bakes the version into the wheel filename, so the
        # docs build writes a manifest rather than us hardcoding it here.
        manifest = json.loads(await (await pyfetch(base + "manifest.json")).string())
        await micropip.install(base + manifest["divtel"])
    return


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import numpy as np
    import astropy.units as u
    import matplotlib.pyplot as plt

    from divtel.telescope import Telescope, Array

    return Array, Telescope, mo, np, plt, u


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # Divergent pointing

    In divergent mode each telescope is offset from the array's mean
    direction, trading photon statistics per shower for a wider combined
    field of view. The **divergence** parameter runs from `0` (all
    telescopes parallel) to `1` (fully radial).
    """)
    return


@app.cell(hide_code=True)
def _():
    return


@app.cell
def _(mo):
    div = mo.ui.slider(
        0, 1, step=0.01, value=0.3, label="divergence", show_value=True,
        full_width=True,
    )
    alt = mo.ui.slider(
        0, 90, step=1, value=45, label="altitude [deg]", show_value=True,
        full_width=True,
    )
    az = mo.ui.slider(
        0, 360, step=1, value=45, label="azimuth [deg]", show_value=True,
        full_width=True,
    )
    mo.vstack([div, alt, az])
    return alt, az, div


@app.cell(hide_code=True)
def _(Array, Telescope, u):
    def hess_1():
        """The four HESS-1 telescopes, on a 100 m square."""
        return Array(
            [
                Telescope(100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
                Telescope(0 * u.m, 100 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
                Telescope(-100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
                Telescope(0 * u.m, -100 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            ]
        )

    array = hess_1()
    return (array,)


@app.cell(hide_code=True)
def _(alt, array, az, div, u):
    # Marimo is reactive: touching a slider re-runs only the cells that read
    # it. This one repoints the array and hands it on under a new name, which
    # is what puts the two plot cells below downstream of the sliders.
    array.divergent_pointing(div.value, alt.value * u.deg, az.value * u.deg)
    pointed = array
    return (pointed,)


@app.cell(hide_code=True)
def _(plt, pointed):
    def _plot(array):
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), layout="constrained")
        for ax, projection in zip(axes, ("xz", "xy", "yz")):
            array.display_2d(projection=projection, ax=ax)
        return fig

    _plot(pointed)
    return


@app.cell(hide_code=True)
def _(np, plt, pointed):
    def _plot3d(array):
        fig = plt.figure(figsize=(5.5, 5))
        ax = fig.add_subplot(111, projection="3d")

        positions = array.positions_array
        pointing_vectors = array.pointing_vectors
        x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
        span = np.ptp(positions, axis=0).max()

        ax.quiver(
            x, y, z,
            pointing_vectors[:, 0],
            pointing_vectors[:, 1],
            pointing_vectors[:, 2],
            length=span,
            color="black",
        )
        ax.scatter(x, y, z, color="tab:blue", s=25)

        ax.set_xlim(x.min() - span, x.max() + span)
        ax.set_ylim(y.min() - span, y.max() + span)
        ax.set_zlim(0, 2 * span)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        return fig

    _plot3d(pointed)
    return


if __name__ == "__main__":
    app.run()

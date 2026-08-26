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
    import marimo as mo
    import numpy as np
    import astropy.units as u
    import matplotlib.pyplot as plt

    from divtel.telescope import Telescope, Array
    from divtel.visualization import display_hyper_fov

    # Render figures as SVG. marimo's PNG path stamps an explicit pixel width
    # on the image -- figsize x 100 -- which overflows a frame narrower than
    # the figure. SVG carries no such width, so it scales to its container.
    plt.rcParams["savefig.format"] = "svg"
    return Array, Telescope, display_hyper_fov, mo, np, plt, u


@app.cell(hide_code=True)
def _(mo):
    # Belt and braces alongside the SVG format above: nothing in the output
    # should be able to out-run its container when embedded in a frame.
    mo.Html(
        """<style>
          img, svg { max-width: 100%; height: auto; }
        </style>"""
    )
    return


@app.cell
def _(mo):
    div = mo.ui.slider(
        0, 1, step=0.005, value=0.02, label="divergence", show_value=True,
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
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.2), layout="constrained")
        for ax, projection in zip(axes, ("xz", "xy", "yz")):
            array.display_2d(projection=projection, ax=ax)
        return fig

    _plot(pointed)
    return


@app.cell(hide_code=True)
def _(display_hyper_fov, np, plt, pointed, u):
    def _plot_sky(array):
        # One figure rather than two cells: the demo is embedded in a frame of
        # fixed height, so the panels have to share a row to stay in view.
        fig = plt.figure(figsize=(11, 4.2), layout="constrained")
        ax3d = fig.add_subplot(1, 2, 1, projection="3d")
        ax_fov = fig.add_subplot(1, 2, 2)

        # Axes are labelled in metres, so hand matplotlib bare numbers.
        positions = array.positions_array.to_value(u.m)
        pointing_vectors = array.pointing_vectors
        x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
        span = np.ptp(positions, axis=0).max()

        ax3d.quiver(
            x, y, z,
            pointing_vectors[:, 0],
            pointing_vectors[:, 1],
            pointing_vectors[:, 2],
            length=span,
            color="black",
        )
        ax3d.scatter(x, y, z, color="tab:blue", s=25)
        ax3d.set_xlim(x.min() - span, x.max() + span)
        ax3d.set_ylim(y.min() - span, y.max() + span)
        ax3d.set_zlim(0, 2 * span)
        ax3d.set_xlabel("x [m]")
        ax3d.set_ylabel("y [m]")
        ax3d.set_zlabel("z [m]")
        ax3d.set_title("pointing on the ground")

        display_hyper_fov(array, ax=ax_fov, m_cut=2)
        return fig

    _plot_sky(pointed)
    return


if __name__ == "__main__":
    app.run()

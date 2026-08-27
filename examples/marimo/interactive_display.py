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
    import marimo as mo
    import numpy as np
    import astropy.units as u
    import matplotlib.pyplot as plt

    from divtel import pointing
    from divtel.telescope import Telescope, Array
    from divtel.visualization import display_hyper_fov

    # SVG scales to its container; marimo's PNG path stamps a fixed pixel
    # width that overflows a frame narrower than the figure.
    plt.rcParams["savefig.format"] = "svg"
    return Array, Telescope, display_hyper_fov, mo, np, plt, pointing, u


@app.cell(hide_code=True)
def _(mo):
    # Belt and braces: caps output size so it can't overflow an embedded frame.
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
def _(alt, array, az, div, pointing, u):
    # Renaming to `pointed` puts the plot cells below downstream of the
    # sliders, since marimo re-runs only cells that read a changed value.
    array.divergent_pointing(div.value, alt.value * u.deg, az.value * u.deg)
    pointed = array
    # G is where the telescopes' pointing directions all trace back to.
    # At div=0 they're parallel and G recedes to infinity, so there's
    # nothing to plot.
    g_point = (
        pointing.pointG_position(
            pointed.barycenter, div.value, alt.value * u.deg, az.value * u.deg
        )
        if div.value > 0
        else None
    )
    return g_point, pointed


@app.cell(hide_code=True)
def _(g_point, mo, np, pointed, u):
    # At low divergence G falls outside the plotted extent below and would
    # seem to have vanished, so print its position and distance regardless.
    if g_point is None:
        _text = "**G**: undefined at zero divergence, pointing is parallel and never converges."
    else:
        _gx, _gy, _gz = g_point.to_value(u.m)
        _dist = np.linalg.norm((g_point - pointed.barycenter).to_value(u.m))
        _text = (
            f"**G** = ({_gx:.0f}, {_gy:.0f}, {_gz:.0f}) m, "
            f"{_dist:.0f} m from the array's barycenter. "
            "Marked with a red x below when it falls within the plotted range."
        )
    mo.md(_text)
    return


@app.cell(hide_code=True)
def _(g_point, plt, pointed, u):
    # Index pairs matching the (x, y) columns Array.display_2d picks for
    # each projection, so G lands on the same axes as the telescopes.
    _axes_for = {"xz": (0, 2), "xy": (1, 0), "yz": (1, 2)}

    def _plot(array, g):
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.2), layout="constrained")
        for ax, projection in zip(axes, ("xz", "xy", "yz")):
            array.display_2d(projection=projection, ax=ax)
            if g is not None:
                # Pin the limits display_2d already chose, so plotting G
                # far out at low divergence can't rescale the array to a
                # dot. G just drops out of frame instead.
                xlim, ylim = ax.get_xlim(), ax.get_ylim()
                i, j = _axes_for[projection]
                gx, gy = g.to_value(u.m)[i], g.to_value(u.m)[j]
                ax.scatter(gx, gy, marker="x", color="red", label="G")
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
        # Marker styles are identical across the three panels, so one legend
        # (on the first) is enough to identify them all.
        axes[0].legend(fontsize="small")
        return fig

    _plot(pointed, g_point)
    return


@app.cell(hide_code=True)
def _(display_hyper_fov, np, plt, pointed, u):
    def _plot_sky(array):
        # One figure, not two cells, so the panels share a row and fit the
        # fixed-height frame this demo is embedded in.
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

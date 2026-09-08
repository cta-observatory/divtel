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
    # Not displayed here -- shown to the right of the three projections
    # below, so dragging a slider and watching its effect stay in the same
    # glance instead of scrolling between two rows.
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
def _(alt, az, div, g_point, mo, plt, pointed, u):
    # Index pairs matching the (x, y) columns Array.display_2d picks for
    # each projection, so G lands on the same axes as the telescopes.
    _axes_for = {"xz": (0, 2), "xy": (1, 0), "yz": (1, 2)}

    def _plot(array, g):
        # Narrower than the row it sits in -- the sliders take the rest of
        # the width, to its right.
        fig, axes = plt.subplots(1, 3, figsize=(7.5, 3.2), layout="constrained")
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

    mo.hstack(
        [_plot(pointed, g_point), mo.vstack([div, alt, az])],
        widths=[2, 1],
        align="center",
    )
    return


@app.cell(hide_code=True)
def _(alt, az, display_hyper_fov, g_point, mo, np, plt, pointed, pointing, u):
    def _plot_3d(array, g, alt_mean, az_mean):
        # Its own figure, not a subplot sharing a row with the hyper-FoV
        # panel: that panel forces an equal-aspect box via
        # `adjustable="datalim"`, so squeezing it into whatever height the
        # dense 3D panel happens to need stretched its altitude axis into a
        # mostly-empty range. Separate figures let each keep its own aspect.
        fig = plt.figure(figsize=(9, 8), layout="constrained")
        ax3d = fig.add_subplot(111, projection="3d")

        positions = array.positions_array.to_value(u.m)
        x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
        barycenter = array.barycenter.to_value(u.m)
        pointing_vectors = array.pointing_vectors
        span = max(np.ptp(positions, axis=0).max(), 50.0)
        arrow_len = span * 0.9
        mean_vec = pointing.alt_az_to_vector(alt_mean, az_mean)

        # Fixed frame: constant limits, sized for this array's 100 m square,
        # so the grid holds still as the sliders move the arrows and G
        # within it instead of resizing around them.
        xy_half = 480.0
        z_lo, z_hi = -750.0, 180.0
        bx, by, bz = barycenter
        lo = np.array([bx - xy_half, by - xy_half, z_lo])
        hi = np.array([bx + xy_half, by + xy_half, z_hi])
        ax3d.set_xlim(lo[0], hi[0])
        ax3d.set_ylim(lo[1], hi[1])
        ax3d.set_zlim(lo[2], hi[2])
        ax3d.set_box_aspect(hi - lo)

        # Ground plane, for a horizon to judge tilt against.
        grid_x, grid_y = np.meshgrid([lo[0], hi[0]], [lo[1], hi[1]])
        grid_z = np.full_like(grid_x, z.mean())
        ax3d.plot_surface(grid_x, grid_y, grid_z, color="0.65", alpha=0.55,
                          linewidth=0, shade=False)

        # Individual telescope pointings.
        ax3d.quiver(x, y, z,
                    pointing_vectors[:, 0], pointing_vectors[:, 1], pointing_vectors[:, 2],
                    length=arrow_len, color="black", linewidth=2.2,
                    arrow_length_ratio=0.15)
        ax3d.scatter(x, y, z, color="tab:blue", s=110, depthshade=False,
                    label="telescopes")
        for xi, yi, zi, tel in zip(x, y, z, array.telescopes):
            ax3d.text(xi, yi, zi, f"  T{tel.id}", fontsize=12, color="tab:blue")

        # Mean pointing, from the barycenter.
        ax3d.quiver(*barycenter, *mean_vec, length=arrow_len * 1.3, color="tab:red",
                    linewidth=3.0, arrow_length_ratio=0.12)
        ax3d.scatter(*barycenter, color="tab:red", marker="+", s=280, linewidths=3.2,
                    label="barycenter / mean pointing")

        # Plumb line through the barycenter, for a vertical to judge tilt against.
        ax3d.plot([bx, bx], [by, by], [z_lo, z_hi], color="steelblue",
                  linewidth=1.8, linestyle=(0, (8, 6)), alpha=0.8)

        # G, and the ray from it through each telescope. Below about
        # div=0.15, G recedes past the bottom of this frame.
        if g is not None:
            gx_, gy_, gz_ = g.to_value(u.m)
            for xi, yi, zi in zip(x, y, z):
                ax3d.plot([gx_, xi], [gy_, yi], [gz_, zi], color="tab:orange",
                         linestyle="--", linewidth=1.6, alpha=0.7)
            ax3d.plot([gx_, bx], [gy_, by], [gz_, bz], color="tab:orange",
                     linestyle="--", linewidth=2.4, alpha=0.9)
            ax3d.scatter([gx_], [gy_], [gz_], color="tab:orange", marker="D", s=150,
                        label="G")

        ax3d.set_xlabel("x [m]", fontsize=13, labelpad=12)
        ax3d.set_ylabel("y [m]", fontsize=13, labelpad=12)
        ax3d.set_zlabel("z [m]", fontsize=13, labelpad=6)
        ax3d.tick_params(labelsize=11)
        ax3d.set_title("divergence geometry", fontsize=16)
        ax3d.legend(fontsize="small", loc="upper left")
        return fig

    def _plot_fov(array):
        fig, ax = plt.subplots(figsize=(6.5, 6.5), layout="constrained")
        display_hyper_fov(array, ax=ax, min_telescopes=2)
        return fig

    mo.hstack(
        [
            _plot_3d(pointed, g_point, alt.value * u.deg, az.value * u.deg),
            _plot_fov(pointed),
        ],
        widths=[1.3, 1],
        align="center",
    )
    return


if __name__ == "__main__":
    app.run()

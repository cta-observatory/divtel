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

    # SVG scales to its container; marimo's PNG path stamps a fixed pixel
    # width that overflows a frame narrower than the figure.
    plt.rcParams["savefig.format"] = "svg"
    return Array, Telescope, mo, np, plt, pointing, u


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
        # Divergence geometry in 3D

        Four telescopes, one mean pointing direction, and the geometry
        `divergent_pointing` builds from them: a virtual point **G** behind
        the array, and one telescope pointing per ray from **G** through
        that telescope's ground position.

        - **black arrows** — each telescope's actual pointing
        - **red arrow** — the mean pointing the array was asked to diverge from
        - **orange diamond + dashed rays** — **G**, and the line from it to
          each telescope; a telescope's pointing is exactly this line's
          direction, continued past the telescope
        """
    )
    return


@app.cell
def _(mo):
    div = mo.ui.slider(
        0, 0.3, step=0.005, value=0.15, label="divergence", show_value=True,
        full_width=True,
    )
    alt = mo.ui.slider(
        0, 90, step=1, value=70, label="mean altitude [deg]", show_value=True,
        full_width=True,
    )
    az = mo.ui.slider(
        0, 360, step=1, value=0, label="mean azimuth [deg]", show_value=True,
        full_width=True,
    )
    mo.vstack([div, alt, az])
    return alt, az, div


@app.cell(hide_code=True)
def _(Array, Telescope, u):
    def four_telescopes():
        """Four telescopes on a 100 m square, HESS-1 style."""
        return Array(
            [
                Telescope(100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
                Telescope(0 * u.m, 100 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
                Telescope(-100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
                Telescope(0 * u.m, -100 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            ]
        )

    array = four_telescopes()
    return (array,)


@app.cell(hide_code=True)
def _(alt, array, az, div, pointing, u):
    alt_mean = alt.value * u.deg
    az_mean = az.value * u.deg
    # Renaming to `pointed` puts the plot cell below downstream of the
    # sliders, since marimo re-runs only cells that read a changed value.
    array.divergent_pointing(div.value, alt_mean, az_mean)
    pointed = array
    # G is where every telescope's pointing traces back to. At div=0 the
    # telescopes point in parallel and G recedes to infinity, so there is
    # nothing to place.
    g_point = (
        pointing.pointG_position(pointed.barycenter, div.value, alt_mean, az_mean)
        if div.value > 0
        else None
    )
    mean_vector = pointing.alt_az_to_vector(alt_mean, az_mean)
    return alt_mean, az_mean, g_point, mean_vector, pointed


@app.cell(hide_code=True)
def _(alt_mean, az_mean, div, g_point, mo, np, pointed, u):
    if g_point is None:
        _text = (
            f"**div = 0** &rarr; **G** is infinitely far below the array; "
            f"all telescopes point straight at (alt={alt_mean:.0f}, az={az_mean:.0f})."
        )
    else:
        _gx, _gy, _gz = g_point.to_value(u.m)
        _norm = 100 * u.m / np.tan(np.arcsin(div.value))
        _dist = np.linalg.norm((g_point - pointed.barycenter).to_value(u.m))
        _text = (
            f"**div = {div.value:.3f}**, mean pointing (alt={alt_mean:.0f}, az={az_mean:.0f})\n\n"
            f"norm = 100 m / tan(asin(div)) = **{_norm:.0f}**\n\n"
            f"**G** = ({_gx:.0f}, {_gy:.0f}, {_gz:.0f}) m, "
            f"**{_dist:.0f} m** from the array's barycenter, behind it along the mean pointing"
        )
    mo.md(_text)
    return


@app.cell(hide_code=True)
def _(g_point, mean_vector, np, plt, pointed, u):
    def _plot(array, g, mean_vec):
        fig = plt.figure(figsize=(8, 7), layout="constrained")
        ax = fig.add_subplot(111, projection="3d")

        positions = array.positions_array.to_value(u.m)
        x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
        barycenter = array.barycenter.to_value(u.m)
        pointing_vectors = array.pointing_vectors
        span = max(np.ptp(positions, axis=0).max(), 50.0)
        arrow_len = span * 0.9

        # The frame is fixed -- constant limits, chosen once for this array
        # and its slider ranges -- rather than refit to wherever the sliders
        # currently put the arrows and G. A refit-every-frame box makes the
        # grid itself appear to slide and rescale as you drag a slider, which
        # reads as the array moving; a fixed box holds still and lets the
        # content move within it instead. The trade-off: below about
        # div=0.15, G recedes past the bottom of frame (the text panel above
        # still reports exactly where it went).
        xy_half = 480.0
        z_lo, z_hi = -750.0, 180.0
        lo = np.array([barycenter[0] - xy_half, barycenter[1] - xy_half, z_lo])
        hi = np.array([barycenter[0] + xy_half, barycenter[1] + xy_half, z_hi])

        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_zlim(lo[2], hi[2])
        ax.set_box_aspect(hi - lo)

        # A ground plane under the array, for a horizon to judge tilt against.
        grid_x, grid_y = np.meshgrid([lo[0], hi[0]], [lo[1], hi[1]])
        grid_z = np.full_like(grid_x, z.mean())
        ax.plot_surface(grid_x, grid_y, grid_z, color="0.65", alpha=0.55,
                        linewidth=0, shade=False)

        # Individual telescope pointings.
        ax.quiver(x, y, z,
                  pointing_vectors[:, 0], pointing_vectors[:, 1], pointing_vectors[:, 2],
                  length=arrow_len, color="black", linewidth=1.6,
                  arrow_length_ratio=0.15)
        ax.scatter(x, y, z, color="tab:blue", s=70, depthshade=False, label="telescopes")
        for xi, yi, zi, tel in zip(x, y, z, array.telescopes):
            ax.text(xi, yi, zi, f"  T{tel.id}", fontsize=8, color="tab:blue")

        # Mean pointing, from the barycenter.
        ax.quiver(*barycenter, *mean_vec, length=arrow_len * 1.3, color="tab:red",
                  linewidth=2.2, arrow_length_ratio=0.12)
        ax.scatter(*barycenter, color="tab:red", marker="+", s=180, linewidths=2.5,
                   label="barycenter / mean pointing")

        # A vertical reference line through the barycenter, top to bottom of
        # frame -- a plumb line to judge tilt against, without the clutter of
        # X/Y arrows that mostly just point into the ground plane anyway.
        bx, by, bz = barycenter
        ax.plot([bx, bx], [by, by], [z_lo, z_hi], color="steelblue",
               linewidth=1.2, linestyle=(0, (8, 6)), alpha=0.8)
        ax.text(bx, by, z_hi, "Z", color="steelblue", fontsize=9,
               fontweight="bold")

        # G, the ray from it through each telescope, and the line back to
        # the barycenter. Below the reference divergence the frame was sized
        # for, G and these lines simply run off the edge.
        if g is not None:
            gx_, gy_, gz_ = g.to_value(u.m)
            for xi, yi, zi in zip(x, y, z):
                ax.plot([gx_, xi], [gy_, yi], [gz_, zi], color="tab:orange",
                       linestyle="--", linewidth=1.1, alpha=0.7)
            ax.plot([gx_, bx], [gy_, by], [gz_, bz], color="tab:orange",
                   linestyle="--", linewidth=1.6, alpha=0.9)
            ax.scatter([gx_], [gy_], [gz_], color="tab:orange", marker="D", s=80,
                      label="G")

        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.legend(fontsize="small", loc="upper left")
        return fig

    _plot(pointed, g_point, mean_vector)
    return


if __name__ == "__main__":
    app.run()

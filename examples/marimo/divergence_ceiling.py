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

    from divtel import pointing as ptg
    from divtel.layout import load_array
    from divtel.visualization import display_hyper_fov

    # SVG scales to its container; marimo's PNG path stamps a fixed pixel
    # width that overflows a frame narrower than the figure.
    plt.rcParams["savefig.format"] = "svg"
    return display_hyper_fov, files, load_array, mo, np, plt, ptg, u


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
    # The divergence ceiling

    Divergent pointing spreads an array's cameras apart so they cover more
    sky. The trade is not free: a shower only one telescope sees cannot be
    reconstructed stereoscopically, so what matters is not the covered area
    but the area seen by two or more telescopes at once, the *stereoscopic
    hyper field of view*. Pushed far enough, that quantity stops rising and
    comes back down — the mean number of telescopes on a typical direction
    falls below two, and further spreading trades stereoscopic sky for
    single-telescope sky the array cannot use. That turning point is the
    *divergence ceiling*, and this notebook locates it for both CTAO sites.

    Each site pairs two telescope types with different cameras: La Palma
    puts four Large-Sized Telescopes inside nine Medium-Sized ones, Paranal
    fourteen Medium-Sized inside thirty-seven Small-Sized. Because the two
    types carry different fields of view and sit at different distances
    from the array centre, they saturate at different divergences, and the
    array as a whole does not simply add their two responses together.
    """)
    return


@app.cell(hide_code=True)
def _(np, u):
    def cap_area_deg2(radius):
        """Exact solid angle of a camera disc of angular radius `radius`, in deg²."""
        steradians = 2 * np.pi * (1 - np.cos(radius.to_value(u.rad)))
        return steradians * (180 / np.pi) ** 2

    def omega(array):
        """Total camera solid angle Ω: the sum of every disc, independent of overlap."""
        return sum(cap_area_deg2(tel.fov_radius) for tel in array.telescopes) * u.deg ** 2

    def parallel_reach(array):
        """Angular radius of the widest camera — what the array sees pointed conventionally."""
        return max(tel.fov_radius for tel in array.telescopes).to(u.deg)

    def angular_spread(array):
        """Smallest and largest angle, in degrees, between a telescope's pointing
        and the array's mean pointing."""
        mean_alt, mean_az = array.mean_pointing
        mean_vector = ptg.alt_az_to_vector(mean_alt, mean_az)
        cos_offset = np.clip(array.pointing_vectors @ mean_vector, -1, 1)
        offsets = np.degrees(np.arccos(cos_offset))
        return offsets.min(), offsets.max()

    def find_ceiling_div(divs, means):
        """The div at which mean multiplicity crosses two, by interpolation."""
        order = np.argsort(means)
        return float(np.interp(2.0, means[order], divs[order]))

    def sweep(array, types, divs, altitudes, az=180 * u.deg):
        """Stereoscopic hyper FoV and mean multiplicity across a div grid.

        Returns two dicts: the whole array at each altitude, keyed by altitude
        in degrees, and each type at the zenith, keyed by type name.
        """
        whole = {}
        for alt in altitudes:
            areas, means = [], []
            for div in divs:
                array.divergent_pointing(float(div), alt, az)
                area, patches = array.hyper_fov(min_telescopes=2)
                mean, _ = array.multiplicity_moments(patches=patches)
                areas.append(area.to_value(u.deg ** 2))
                means.append(mean)
            whole[int(round(alt.to_value(u.deg)))] = (np.array(areas), np.array(means))

        zenith = 90 * u.deg
        by_type = {}
        for name, sub in array.group_by(types).items():
            areas, means = [], []
            for div in divs:
                sub.divergent_pointing(float(div), zenith, az)
                area, patches = sub.hyper_fov(min_telescopes=2)
                mean, _ = sub.multiplicity_moments(patches=patches)
                areas.append(area.to_value(u.deg ** 2))
                means.append(mean)
            by_type[name] = (np.array(areas), np.array(means))

        return whole, by_type

    def ceiling_report(array, alt, az=180 * u.deg, lo=0.001, hi=0.15, n=80):
        """Everything Table 1 needs, computed at one reference altitude."""
        reach = parallel_reach(array)
        area1 = cap_area_deg2(reach) * u.deg ** 2
        total_omega = omega(array)
        mean_parallel = float(total_omega / area1)

        divs = np.linspace(lo, hi, n)
        means = []
        for div in divs:
            array.divergent_pointing(float(div), alt, az)
            _, patches = array.hyper_fov(min_telescopes=1)
            mean, _ = array.multiplicity_moments(patches=patches)
            means.append(mean)
        ceiling_div = find_ceiling_div(divs, np.array(means))

        array.divergent_pointing(ceiling_div, alt, az)
        area2, patches = array.hyper_fov(min_telescopes=2)
        _, spread = angular_spread(array)
        multiplicity, per_multiplicity = array.multiplicity_profile(patches=patches)
        total = per_multiplicity.sum()
        waste = float(per_multiplicity[multiplicity == 1].sum() / total) if 1 in multiplicity else 0.0
        deep = float(per_multiplicity[multiplicity >= 3].sum() / total)

        return {
            "reach": reach.to_value(u.deg),
            "area1": area1.to_value(u.deg ** 2),
            "omega": total_omega.to_value(u.deg ** 2),
            "mean_parallel": mean_parallel,
            "ceiling_div": ceiling_div,
            "spread": spread,
            "area2": area2.to_value(u.deg ** 2),
            "gain": float(area2 / area1),
            "omega_frac": float(area2 / (total_omega / 2)),
            "waste": waste,
            "deep": deep,
        }

    return angular_spread, cap_area_deg2, ceiling_report, find_ceiling_div, omega, parallel_reach, sweep


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## CTAO-North: La Palma

    Four LSTs sit inside nine MSTs. An LST camera has a field-of-view
    radius of about 2.15°, an MST's about 3.84°, so the outer, more
    numerous subarray also carries the wider camera.
    """)
    return


@app.cell(hide_code=True)
def _(files, load_array):
    ARRAY_N = load_array(
        files("divtel") / "data" / "cta-north-lapalma-alpha-prod6.ecsv"
    )
    # Ids follow the CTAO convention: LSTs 1-4, MSTs 5-14.
    TYPES_N = {"LST": range(1, 5), "MST": range(5, 14)}
    return ARRAY_N, TYPES_N


@app.cell(hide_code=True)
def _(ARRAY_N, TYPES_N, np, sweep, u):
    divs_N = np.linspace(0.002, 0.16, 40)
    alts_N = np.array([30, 50, 70, 90]) * u.deg
    whole_N, by_type_N = sweep(ARRAY_N, TYPES_N, divs_N, alts_N)
    return alts_N, by_type_N, divs_N, whole_N


@app.cell(hide_code=True)
def _(by_type_N, divs_N, plt, whole_N):
    def _fig3():
        fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.5))
        for alt, (areas, means) in whole_N.items():
            line, = left.plot(divs_N, areas, label=f"{alt}° (whole array)")
            right.plot(divs_N, means, color=line.get_color())
        for name, (areas, means) in by_type_N.items():
            line, = left.plot(divs_N, areas, "--", label=f"{name} (zenith)")
            right.plot(divs_N, means, "--", color=line.get_color())
        right.axhline(2, color="k", linewidth=1, linestyle=":", label="stereoscopic floor")

        left.set_xlabel("div")
        left.set_ylabel("stereoscopic hyper FoV [deg$^2$]")
        right.set_xlabel("div")
        right.set_ylabel("mean multiplicity")
        left.grid(alpha=0.3)
        right.grid(alpha=0.3)
        left.legend(frameon=False, fontsize=8)
        right.legend(frameon=False, fontsize=8)
        fig.suptitle("CTAO-North: stereoscopic coverage and mean multiplicity against div")
        fig.tight_layout()
        return fig

    _fig3()
    return


@app.cell(hide_code=True)
def _(by_type_N, divs_N, mo, np, whole_N):
    def _prose():
        peak_area = {alt: areas.max() for alt, (areas, _) in whole_N.items()}
        peak_div = {alt: divs_N[np.argmax(areas)] for alt, (areas, _) in whole_N.items()}
        span_lo, span_hi = min(peak_area.values()), max(peak_area.values())
        lo_alt, hi_alt = min(whole_N), max(whole_N)

        lst_areas, _ = by_type_N["LST"]
        mst_areas, _ = by_type_N["MST"]
        lst_div, lst_peak = divs_N[np.argmax(lst_areas)], lst_areas.max()
        mst_div, mst_peak = divs_N[np.argmax(mst_areas)], mst_areas.max()

        return mo.md(f"""
        Each type saturates on its own schedule. The four LSTs, the
        narrower camera, peak near div = {lst_div:.3f} at {lst_peak:.0f}
        deg²; the nine MSTs peak later, near div = {mst_div:.3f}, at
        {mst_peak:.0f} deg². The whole-array curve is not their sum: at low
        div it tracks the MSTs, then pulls ahead of them once the MSTs stop
        overlapping each other and cross-type pairs take over.

        Altitude moves the peak without moving its height by much. Between
        {lo_alt}° and {hi_alt}° altitude the div that maximises coverage
        shifts from {peak_div[hi_alt]:.3f} to {peak_div[lo_alt]:.3f}, while
        the coverage it buys stays within
        {(span_hi - span_lo) / span_hi * 100:.0f}% of {span_hi:.0f} deg².
        The parameter has to be chosen for the pointing; what it delivers
        is close to a property of the array. On the right, mean
        multiplicity falls monotonically and crosses the stereoscopic
        floor of two at a div the next section pins down exactly.
        """)

    _prose()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## CTAO-South: Paranal

    Fourteen MSTs sit inside thirty-seven SSTs. The SST mirror is the
    smaller of the two, but its shorter focal length gives it the wider
    field of view: about 4.40° against 3.75° for the MST. Same shape as
    La Palma — the outer, more numerous subarray carries the wider camera
    — but the two cameras here are much closer in size.
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
def _(ARRAY_S, TYPES_S, np, sweep, u):
    divs_S = np.linspace(0.002, 0.16, 40)
    alts_S = np.array([30, 50, 70, 90]) * u.deg
    whole_S, by_type_S = sweep(ARRAY_S, TYPES_S, divs_S, alts_S)
    return alts_S, by_type_S, divs_S, whole_S


@app.cell(hide_code=True)
def _(by_type_S, divs_S, plt, whole_S):
    def _fig3():
        fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.5))
        for alt, (areas, means) in whole_S.items():
            line, = left.plot(divs_S, areas, label=f"{alt}° (whole array)")
            right.plot(divs_S, means, color=line.get_color())
        for name, (areas, means) in by_type_S.items():
            line, = left.plot(divs_S, areas, "--", label=f"{name} (zenith)")
            right.plot(divs_S, means, "--", color=line.get_color())
        right.axhline(2, color="k", linewidth=1, linestyle=":", label="stereoscopic floor")

        left.set_xlabel("div")
        left.set_ylabel("stereoscopic hyper FoV [deg$^2$]")
        right.set_xlabel("div")
        right.set_ylabel("mean multiplicity")
        left.grid(alpha=0.3)
        right.grid(alpha=0.3)
        left.legend(frameon=False, fontsize=8)
        right.legend(frameon=False, fontsize=8)
        fig.suptitle("CTAO-South: stereoscopic coverage and mean multiplicity against div")
        fig.tight_layout()
        return fig

    _fig3()
    return


@app.cell(hide_code=True)
def _(by_type_S, divs_S, mo, np, whole_S):
    def _prose():
        peak_area = {alt: areas.max() for alt, (areas, _) in whole_S.items()}
        peak_div = {alt: divs_S[np.argmax(areas)] for alt, (areas, _) in whole_S.items()}
        span_lo, span_hi = min(peak_area.values()), max(peak_area.values())
        lo_alt, hi_alt = min(whole_S), max(whole_S)

        mst_areas, _ = by_type_S["MST"]
        sst_areas, _ = by_type_S["SST"]
        mst_div, mst_peak = divs_S[np.argmax(mst_areas)], mst_areas.max()
        sst_div, sst_peak = divs_S[np.argmax(sst_areas)], sst_areas.max()

        return mo.md(f"""
        Same rise-then-collapse shape as La Palma, just bigger and faster.
        The SSTs, more numerous and wider-eyed, peak near div = {sst_div:.3f}
        at {sst_peak:.0f} deg²; the MSTs peak lower, near div = {mst_div:.3f},
        at {mst_peak:.0f} deg². At low div the whole array tracks the SSTs
        almost exactly — the MSTs' narrower sky sits inside the SSTs' wider
        one, so folding them in adds multiplicity rather than area — until
        cross-type pairs start carrying coverage the SSTs alone have lost.

        Across {lo_alt}° to {hi_alt}° altitude the div that maximises
        coverage shifts from {peak_div[hi_alt]:.3f} to
        {peak_div[lo_alt]:.3f}, while the coverage it buys stays within
        {(span_hi - span_lo) / span_hi * 100:.0f}% of {span_hi:.0f} deg².
        With fifty-one telescopes of two similarly sized cameras to pair
        up, Paranal keeps mean multiplicity above the stereoscopic floor
        of two well past the div where either subarray alone has given
        out.
        """)

    _prose()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Locating the ceiling

    Two numbers answer how far an array can usefully be spread. The
    *parallel reach* is the angular radius of its widest camera — what the
    array sees pointed conventionally, however many telescopes it has. The
    *ceiling divergence* is the div at which mean multiplicity falls to
    two. Below it, spreading converts sky seen many times into sky seen
    twice, and the stereoscopic footprint $A(2)$ grows; above it the
    typical direction is seen by fewer than two telescopes, and spreading
    converts stereoscopic sky into single-telescope sky the array cannot
    reconstruct from. $A(2)$ turns over exactly at the ceiling and falls
    back towards zero past it — asked to spread further, the array
    complies and returns a configuration whose stereoscopic coverage is
    smaller than at the ceiling, not larger.

    Divergence conserves the total camera solid angle $\Omega$, the sum of
    every disc regardless of overlap, so an array able to place exactly
    two telescopes on every direction it covers would reach
    $A(2) = \Omega / 2$. That bound is never attained, because the
    multiplicity distribution at the ceiling has spread around its mean of
    two: some sky is seen once and wasted, some three times or more and
    deeper than stereoscopy requires. Both figures below are computed at
    60° altitude, so the two arrays can be read off the same table.
    """)
    return


@app.cell(hide_code=True)
def _(ARRAY_N, ARRAY_S, ceiling_report, u):
    report_N = ceiling_report(ARRAY_N, 60 * u.deg)
    report_S = ceiling_report(ARRAY_S, 60 * u.deg)
    return report_N, report_S


@app.cell(hide_code=True)
def _(mo, report_N, report_S):
    mo.md(f"""
    | | CTAO-North | CTAO-South |
    |---|---|---|
    | Telescopes | 4 LST, 9 MST | 14 MST, 37 SST |
    | Parallel reach | {report_N['reach']:.1f}° | {report_S['reach']:.1f}° |
    | Parallel field of view $A(1)$ | {report_N['area1']:.0f} deg² | {report_S['area1']:.0f} deg² |
    | Parallel mean multiplicity | {report_N['mean_parallel']:.1f} | {report_S['mean_parallel']:.1f} |
    | Total camera solid angle $\\Omega$ | {report_N['omega']:.0f} deg² | {report_S['omega']:.0f} deg² |
    | Ceiling divergence | {report_N['ceiling_div']:.3f} | {report_S['ceiling_div']:.3f} |
    | Spread at the ceiling | {report_N['spread']:.1f}° | {report_S['spread']:.1f}° |
    | Stereoscopic footprint $A(2)$ | {report_N['area2']:.0f} deg² | {report_S['area2']:.0f} deg² |
    | Gain over parallel | {report_N['gain']:.1f}× | {report_S['gain']:.1f}× |
    | Fraction of the $\\Omega/2$ bound | {report_N['omega_frac']:.0%} | {report_S['omega_frac']:.0%} |
    | Covered sky wasted at multiplicity 1 | {report_N['waste']:.0%} | {report_S['waste']:.0%} |
    | Covered sky deeper than needed ($\\ge$3) | {report_N['deep']:.0%} | {report_S['deep']:.0%} |

    *Spread is the largest angle between a telescope's pointing and the
    array's mean pointing, at the ceiling. Camera radii come from the
    prod6 layout files, at 60° altitude, 180° azimuth.*
    """)
    return


@app.cell(hide_code=True)
def _(ARRAY_N, ARRAY_S, display_hyper_fov, plt, report_N, report_S, u):
    def _ceiling_maps():
        fig, (left, right) = plt.subplots(1, 2, figsize=(11, 5))

        ARRAY_N.divergent_pointing(report_N["ceiling_div"], 60 * u.deg, 180 * u.deg)
        display_hyper_fov(ARRAY_N, ax=left, min_telescopes=1, show_area=False)
        left.set_title(f"CTAO-North, div = {report_N['ceiling_div']:.3f}")

        ARRAY_S.divergent_pointing(report_S["ceiling_div"], 60 * u.deg, 180 * u.deg)
        display_hyper_fov(ARRAY_S, ax=right, min_telescopes=1, show_area=False)
        right.set_title(f"CTAO-South, div = {report_S['ceiling_div']:.3f}")

        fig.suptitle("Coverage at the ceiling divergence, shaded by multiplicity")
        fig.tight_layout()
        return fig

    _ceiling_maps()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Two things the table and the maps show

    **The two arrays are alike in reach and differ in depth.** CTAO-South
    is four times larger on the ground with four times the telescopes, yet
    the two parallel reaches sit within half a degree of each other and
    the two ceiling divergences within 0.002 — both quantities depend on
    the ratio of camera size to array scale, and that ratio happens to be
    similar at La Palma and Paranal. What the ceiling buys differs by a
    much larger factor, and the parallel mean multiplicity says why: since
    the gain available over parallel pointing is half of it, South, with
    roughly four times North's parallel depth, converts that depth into a
    correspondingly larger stereoscopic footprint. Divergence does not
    make a large array reach further than a small one; it lets a deep
    array trade depth it does not need for sky, and a shallow array has
    less to trade.

    **Nearly half of what a divergent array watches is wasted.** Neither
    array attains the $\Omega/2$ bound, and the reason is not the mean
    multiplicity but its spread: mean multiplicity is two at the ceiling
    in both cases, but a large share of the covered sky sits at
    multiplicity one, seen by a single telescope and useless for
    stereoscopic reconstruction, while a comparable share sits at three or
    more, deeper than stereoscopy requires. The maps above draw that
    structure directly — an over-covered core where several cameras pile
    up, ringed by a broad band seen by exactly one telescope. A symmetric
    fan is shaped like the array's ground footprint and has no way to move
    depth from where there is too much of it to where there is none.
    Recovering that waste means giving up the single global parameter and
    letting each telescope's pointing follow the target instead of the
    ground.
    """)
    return


if __name__ == "__main__":
    app.run()

#!/usr/bin/env python
"""
Everything the GW170817 study page shows, computed from scratch.

The page quotes a lot of numbers and it must not be possible for any of them to
drift from the code. So nothing on it is typed by hand: this script runs the
strategies, writes the figures, and writes the numbers as reStructuredText
substitutions and table fragments that the page includes. Change the geometry
and the prose changes with it, or the build fails on a substitution that no
longer exists.

It reads the credible regions bundled in ``divtel/data/gw170817``, so it needs
neither the network nor the ``divtel[skymap]`` extra, and ``docs/conf.py`` runs
it before Sphinx reads a source file. Run it alone to see the numbers::

    python docs/scripts/make_gw170817.py --outdir /tmp/gw

Set ``DIVTEL_DOCS_SKIP_GW170817=1`` to skip it during a docs build, which is
worth doing while editing prose and is not worth doing before publishing.
"""

from __future__ import annotations

import argparse
import json
from importlib.resources import files
from pathlib import Path

import astropy.units as u
import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from astropy.coordinates import SkyCoord  # noqa: E402
from astropy.time import Time  # noqa: E402

from divtel import strategy, visualization  # noqa: E402
from divtel.layout import load_array  # noqa: E402
from divtel.observation import altaz_track, observable_windows  # noqa: E402
from divtel.pointing import (  # noqa: E402
    best_pointing,
    div_for_multiplicity,
    pointing_spread,
)
from divtel.region import SkyRegion  # noqa: E402

DATA = files("divtel") / "data"
REGIONS = DATA / "gw170817"

# GW170817 merged at 12:41:04 UTC on 17 August 2017, and the optical counterpart
# SSS17a was found in NGC 4993 at this position.
MERGER = "2017-08-17T12:41:04"
SOURCE = SkyCoord(ra=197.450374 * u.deg, dec=-23.381495 * u.deg)

SITES = {"CTAO-North (La Palma)": "north", "CTAO-South (Paranal)": "south"}
LAYOUTS = {"north": "cta-north-lapalma-alpha-prod6.ecsv",
           "south": "cta-south-paranal-alpha-prod6.ecsv"}

# The alert sequence, and when each circular went out relative to the merger.
ALERTS = ["bayestar_hl", "bayestar_hlv", "lalinference_prelim"]

# Every region is placed at the same elevation, so what is compared is
# localizations and not the elevations they happened to fall at. An array
# looking near the horizon is foreshortened, which changes the spread a given
# divergence produces and the multiplicity it leaves behind.
ALT, AZ = 60 * u.deg, 180 * u.deg
LEVEL, M_CUT = 0.9, 2

# The case every strategy is run against: the first alert map, from the site
# that could actually see it. It is the one the divergence-only study found no
# divergence could cover.
SHOWCASE, SHOWCASE_SITE = "bayestar_hl", "south"

SUBARRAY_COUNTS = [2, 3, 4, 6, 8, 10, 14, 17]
BETAS = [0.0, 0.5, 1.0, 1.5, 2.0]

plt.rcParams.update({"font.size": 9, "figure.dpi": 130,
                     "savefig.bbox": "tight", "axes.spines.top": False,
                     "axes.spines.right": False})


# -- helpers ---------------------------------------------------------------

def load_region(name, level=LEVEL):
    """One bundled credible region, still where the sky put it."""
    return SkyRegion.from_table(REGIONS / f"{name}_{round(100 * level)}.ecsv.gz")


def placed(name, level=LEVEL, alt=ALT, az=AZ):
    """The same region set down over an array, and cut at the horizon."""
    return load_region(name, level).place(alt, az).visible_part(0 * u.deg)


def placed_bands(name, levels=(0.5, 0.9, 0.95), alt=ALT, az=AZ):
    """
    Every credible band of one map, placed together and cut at the horizon.

    One rotation for all of them, taken from the level the study is scored on.
    Each band has its own centroid, so rotating each by its own would slide
    them apart instead of nesting them.
    """
    rotation = load_region(name, LEVEL).placement_rotation(alt, az)
    return {level: load_region(name, level).rotate(rotation)
                                           .visible_part(0 * u.deg)[0]
            for level in levels}


def array(site=SHOWCASE_SITE):
    return load_array(DATA / LAYOUTS[site])


def hours(quantity):
    return float(u.Quantity(quantity).to_value(u.hour))


def since_merger(time):
    return hours((Time(time) - Time(MERGER)).to(u.hour))


def score(pointed, region, exposures=1):
    """The row every strategy is reported as."""
    described = strategy.describe(pointed, region, M_CUT)
    cos_radii = np.cos([t.fov_radius.to_value(u.rad) for t in pointed.telescopes])
    on_target = int(((region.directions @ pointed.pointing_vectors.T) >= cos_radii)
                    .any(axis=0).sum())
    return {
        "covered": described["covered_stereo"],
        "mean": described["mean_covered"],
        "mismatch": described["mismatch"],
        "on_target": on_target,
        "telescopes": len(pointed.telescopes),
        "exposures": exposures,
        "profile": strategy.multiplicity_by_probability(pointed, region),
    }


# -- the study ------------------------------------------------------------

def visibility(args):
    """When, and from where, the localization could be observed at all."""
    # The probability-weighted centre of the map observers actually had in hand
    # when either site got dark, which is where an array would have pointed.
    target = _as_coord(load_region("bayestar_hlv").centroid)

    result = {"target": {"ra_deg": float(target.ra.deg),
                         "dec_deg": float(target.dec.deg)},
              "sites": {}}

    # The whole five days at once, so the reader sees the source rise and set
    # rather than a row of thin spikes. The windows are only an hour or three
    # each out of a hundred and twenty.
    step = 10 * u.min
    times = Time(MERGER) + np.arange(
        int((5 * u.day / step).to_value(u.dimensionless_unscaled))) * step
    elapsed = (times - Time(MERGER)).to_value(u.hour)

    figure, axes = plt.subplots(2, 1, figsize=(9.6, 4.6), sharex=True,
                                sharey=True)
    for ax, (label, site) in zip(axes, SITES.items(), strict=True):
        windows = observable_windows(target, site, MERGER, 5 * u.day,
                                     step=5 * u.min)
        result["sites"][label] = {
            "windows": [{"start": w.start.isot, "hours": hours(w.duration),
                         "alt_max_deg": float(w.alt_max.to_value(u.deg)),
                         "since_merger_h": since_merger(w.start)}
                        for w in windows],
            "total_hours": sum(hours(w.duration) for w in windows),
            "alt_max_deg": max((float(w.alt_max.to_value(u.deg))
                                for w in windows), default=0.0),
            "first_h": since_merger(windows[0].start) if windows else None,
        }

        altitude, sun = altaz_track(target, site, times)
        altitude = altitude.to_value(u.deg)
        dark = sun.to_value(u.deg) <= -18

        ax.fill_between(elapsed, 0, 90, where=dark, color="0.92", linewidth=0,
                        step="mid")
        ax.plot(elapsed, altitude, color="0.65", linewidth=1)
        observable = dark & (altitude >= 0)
        ax.fill_between(elapsed, 0, altitude, where=observable, color="#2a78d6",
                        alpha=0.75, linewidth=0, step="mid")

        ax.axhline(0, color="0.4", linewidth=1)
        ax.set_ylim(-10, 90)
        ax.set_xlim(0, 120)
        ax.set_ylabel("altitude [deg]")
        ax.set_title(
            f"{label} — {result['sites'][label]['total_hours']:.0f} h "
            f"observable, never above "
            f"{result['sites'][label]['alt_max_deg']:.0f}$\\degree$"
            if result["sites"][label]["alt_max_deg"] < 20 else
            f"{label} — {result['sites'][label]['total_hours']:.0f} h "
            f"observable, first window at "
            f"+{result['sites'][label]['first_h']:.1f} h",
            loc="left")

    axes[-1].set_xlabel("hours after the merger")
    axes[0].annotate("dark sky", xy=(0.5, 0.86), xycoords="axes fraction",
                     color="0.45", fontsize=8)
    figure.tight_layout()
    figure.savefig(args.figures / "visibility.png")
    plt.close(figure)
    return result


def _location(site):
    from divtel.observation import _resolve_site

    return _resolve_site(site)


def _as_coord(vector):
    """A unit vector in a celestial frame as an ICRS `SkyCoord`."""
    return SkyCoord(
        ra=np.arctan2(vector[1], vector[0]) * u.rad,
        dec=np.arcsin(np.clip(vector[2], -1.0, 1.0)) * u.rad, frame="icrs")


def maps(args):
    """How big each alert map was, and what that costs an array."""
    result = {}
    levels = [0.5, 0.9, 0.95]
    shades = ["#2a78d6", "#7aacea", "#bfd6f2"]

    # Two rows. The Mollweide says where on the sky each map is, and shows the
    # first one's two lobes a hundred degrees apart; the equal-area close-ups
    # below share a scale, so folding Virgo in is visible as an area rather than
    # only as a number.
    figure = plt.figure(figsize=(4.6 * len(ALERTS), 5.6))
    centre_ra = load_region(ALERTS[0], 0.9).coord.ra.mean()
    all_bands = {}

    for index, name in enumerate(ALERTS):
        bands = {level: load_region(name, level) for level in levels}
        all_bands[name] = bands
        region = bands[LEVEL]
        extent = region.describe_extent()
        result[name] = {
            "label": region.meta["label"],
            "gcn": region.meta["gcn"],
            "since_merger_h": since_merger(region.meta["t_available"]),
            "area_deg2": float(extent["area"].to_value(u.deg**2)),
            "equivalent_radius_deg": float(
                extent["equivalent_radius"].to_value(u.deg)),
            "containment_radius_deg": float(
                extent["containment_radius"].to_value(u.deg)),
            "max_separation_deg": float(extent["max_separation"].to_value(u.deg)),
            "elongation": extent["elongation"],
        }

        ax = figure.add_subplot(2, len(ALERTS), index + 1, projection="mollweide")
        visualization.sky_bands(
            bands, source=SOURCE, centre_ra=centre_ra, ax=ax,
            title=f"{region.meta['label']}, available "
                  f"+{result[name]['since_merger_h']:.1f} h\n"
                  f"{result[name]['area_deg2']:.0f} deg$^2$, "
                  f"{result[name]['max_separation_deg']:.0f}$\\degree$ long")

    # The close-ups, all at the same scale, so the three read against each other.
    half = 1.1 * max(np.abs(visualization.project_region(bands[0.95])).max()
                     for bands in all_bands.values())
    for index, name in enumerate(ALERTS):
        ax = figure.add_subplot(2, len(ALERTS), len(ALERTS) + index + 1)
        for level, shade in zip(levels, shades, strict=True):
            x, y = visualization.project_region(all_bands[name][level])
            ax.scatter(x, y, s=1.6, marker=".", linewidths=0, color=shade)
        ax.set_xlim(-half, half)
        ax.set_ylim(-half / 4, half / 4)
        ax.set_aspect("equal")
        ax.set_xlabel("equal-area, long axis along the page [deg]")

    handles = [plt.Line2D([], [], marker="o", linestyle="none", markersize=6,
                          color=shade, label=f"{level:.0%} credible")
               for level, shade in zip(levels, shades, strict=True)]
    handles.append(plt.Line2D([], [], marker="*", linestyle="none", markersize=10,
                              color="k", label="SSS17a in NGC 4993"))
    figure.tight_layout()
    figure.legend(handles=handles, frameon=False, ncol=4, loc="lower center",
                  bbox_to_anchor=(0.5, -0.03))
    figure.savefig(args.figures / "maps.png")
    plt.close(figure)
    return result


def arrays():
    """What each array reaches pointed conventionally, and where its ceiling is."""
    result = {}
    for site in ("north", "south"):
        telescopes = array(site)
        radii = np.array([t.fov_radius.to_value(u.rad)
                          for t in telescopes.telescopes])
        omega = float((2 * np.pi * (1 - np.cos(radii))).sum()) * u.sr

        ceiling = div_for_multiplicity(telescopes, ALT, AZ, target=2.0)
        telescopes.divergent_pointing(ceiling, ALT, AZ)
        stereo, _ = telescopes.hyper_fov(min_telescopes=2)

        telescopes.divergent_pointing(0.0, ALT, AZ)
        parallel, _ = telescopes.hyper_fov(min_telescopes=2)

        result[site] = {
            "telescopes": len(telescopes.telescopes),
            "parallel_reach_deg": float(
                strategy.stereo_radius(telescopes, M_CUT).to_value(u.deg)),
            "parallel_area_deg2": float(parallel.to_value(u.deg**2)),
            "ceiling_div": ceiling,
            "ceiling_spread_deg": float(
                pointing_spread(telescopes)["max"].to_value(u.deg)),
            "ceiling_stereo_deg2": float(stereo.to_value(u.deg**2)),
            "omega_deg2": float(omega.to(u.deg**2).value),
            "half_omega_deg2": float(omega.to(u.deg**2).value) / 2,
        }
        result[site]["ceiling_fraction"] = (result[site]["ceiling_stereo_deg2"]
                                            / result[site]["half_omega_deg2"])
    return result


def divergence(region, site=SHOWCASE_SITE, steps=13, div_max=0.30):
    """The most any single divergence achieves, aim searched at each."""
    telescopes = array(site)
    rows, best, previous = [], None, None
    for div in np.linspace(0.0, div_max, steps):
        aim = best_pointing(telescopes, region, float(div), M_CUT, start=previous)
        previous = (aim["alt"], aim["az"])
        row = score(telescopes, region)
        row["div"] = float(div)
        row["spread_deg"] = float(
            pointing_spread(telescopes)["max"].to_value(u.deg))
        rows.append(row)
        if best is None or row["covered"] > best["covered"]:
            best = row
    return rows, best


def fewest_groups(rows, slack=0.005):
    """
    The smallest split that gives up nothing worth having.

    More groups is not better. Coverage climbs with the number of groups and
    then flattens, while the telescopes on any given shower keep falling, so
    taking the *best* coverage picks a needlessly thin array. This takes the
    fewest groups that come within `slack` of the best coverage on offer, which
    is the trade an observer would make.
    """
    best = max(row["covered"] for row in rows)
    return min((row for row in rows if row["covered"] >= best - slack),
               key=lambda row: row["count"])


def strategies(args, region, bands):
    """Every way of covering the showcase map, scored the same way."""
    result = {}

    # -- parallel, and the divergence scan ---------------------------------
    scan, best_div = divergence(region)
    result["scan"] = scan
    result["parallel"] = scan[0]
    result["divergent"] = best_div

    # -- sequential tiling -------------------------------------------------
    tiles = strategy.tile_region(array(), region, pointings=10, m_cut=M_CUT)
    result["tiles"] = [{"alt_deg": float(t["alt"].to_value(u.deg)),
                        "az_deg": float(t["az"].to_value(u.deg)),
                        "new": t["new"], "cumulative": t["cumulative"]}
                       for t in tiles]
    # How many ordinary pointings it takes to match one divergent one, and to
    # get essentially all of it.
    result["tiles_to_match_divergence"] = next(
        (n for n, t in enumerate(tiles, 1)
         if t["cumulative"] >= best_div["covered"]), None)
    result["tiles_to_99"] = next(
        (n for n, t in enumerate(tiles, 1) if t["cumulative"] >= 0.99), None)
    telescopes = array()
    telescopes.divergent_pointing(0.0, tiles[0]["alt"], tiles[0]["az"])
    result["tiling"] = score(telescopes, region,
                             exposures=result["tiles_to_match_divergence"])
    result["tiling"]["covered"] = tiles[
        result["tiles_to_match_divergence"] - 1]["cumulative"]

    # -- sub-arrays --------------------------------------------------------
    rows = []
    for count in SUBARRAY_COUNTS:
        telescopes = array()
        summary = strategy.point_subarrays(telescopes, region, count, M_CUT)
        row = score(telescopes, region)
        row["count"] = count
        row["sizes"] = summary["sizes"]
        rows.append(row)
    result["subarray_scan"] = rows
    result["subarrays"] = fewest_groups(rows)

    telescopes = array()
    strategy.point_subarrays(telescopes, region, result["subarrays"]["count"], M_CUT)
    subarray_array = telescopes

    # -- unequal sub-arrays ------------------------------------------------
    rows = []
    for count in SUBARRAY_COUNTS:
        telescopes = array()
        summary = strategy.weighted_split(telescopes, region, count, M_CUT)
        row = score(telescopes, region)
        row["count"] = count
        row["sizes"] = summary["sizes"]
        rows.append(row)
    result["weighted_scan"] = rows
    result["weighted"] = fewest_groups(rows)

    # -- shaped pointing ---------------------------------------------------
    rows = {}
    shaped_array = None
    for beta in BETAS:
        telescopes = array()
        solved = strategy.shaped_pointing(telescopes, region, beta=beta,
                                          m_cut=M_CUT)
        row = score(telescopes, region)
        row["beta"] = beta
        row["efficiency"] = solved["efficiency"]
        row["reach"] = solved["reach"]
        rows[f"{beta}"] = row
        if beta == 1.0:
            shaped_array = telescopes
    result["beta_scan"] = rows
    result["shaped"] = rows["1.0"]
    result["camera_blur"] = strategy.camera_blur(array(), region)

    # -- the pictures ------------------------------------------------------
    panels(args, region, bands, subarray_array, shaped_array, best_div)
    return result


def panels(args, region, bands, subarray_array, shaped_array, best_div):
    """The four strategies over the same region: cameras, depth, and profile."""
    # Both aimed by the same search the scan used, so the panels show the same
    # configurations the table scores. Aiming a parallel array at the region's
    # centroid instead would be worse than it looks: for an arc the weighted
    # centroid lies off the arc, so the array would stare between the lobes.
    parallel = array()
    best_pointing(parallel, region, 0.0, M_CUT)

    diverged = array()
    best_pointing(diverged, region, best_div["div"], M_CUT)

    shown = {
        "parallel": parallel,
        f"divergent, div = {best_div['div']:.3f}": diverged,
        "sub-arrays": subarray_array,
        "shaped, $\\beta$ = 1": shaped_array,
    }

    # One scale down the column, or each panel is normalised to its own maximum
    # and they cannot be compared. The parallel array piles all fifty-one
    # telescopes on one spot, and scaling to that leaves every other panel a
    # uniform pale wash -- so the scale comes from the other three, and the
    # parallel panel is allowed to clip.
    vmax = int(np.percentile(
        np.concatenate([region.multiplicity(pointed)
                        for label, pointed in shown.items()
                        if label != "parallel"]), 99))

    span = visualization.region_span(region)
    figure, axes = plt.subplots(len(shown), 1,
                                figsize=(8.4, 2.2 + 8.4 * span * len(shown)))
    dots = None
    for ax, (label, pointed) in zip(axes, shown.items(), strict=True):
        dots = visualization.multiplicity_over_region(
            pointed, region, ax=ax, title=label, vmax=vmax, colorbar=False,
            size=6.0)
        ax.set_xlabel("")
    axes[-1].set_xlabel("along the region [deg]")
    bar = figure.colorbar(dots, ax=list(axes), fraction=0.02, pad=0.02)
    bar.set_label("telescopes seeing this direction")
    figure.savefig(args.figures / "strategies.png")
    plt.close(figure)

    figure, ax = plt.subplots(figsize=(7.0, 3.6))
    visualization.multiplicity_by_probability(
        {label: strategy.multiplicity_by_probability(pointed, region)
         for label, pointed in shown.items()}, ax=ax, m_cut=M_CUT)
    figure.savefig(args.figures / "gradient.png")
    plt.close(figure)

    footprints(args, region, bands, shown)


def footprints(args, region, bands, shown):
    """
    Where the cameras actually are, over the probability they are covering.

    The shading in ``strategies.png`` is the result; this is the cause. It is
    also the only picture that shows the telescope *types* being sent to
    different parts of the map, which is a degree of freedom neither divergence
    nor an even sub-array split can use.
    """
    frame = visualization.projection_frame(region)

    figure, axes = plt.subplots(len(shown), 1, sharex=True, sharey=True)

    # One scale down the column, framed on the widest configuration, or the
    # panels cannot be compared -- and a divergent array flung across the sky
    # would be cropped to look as tidy as a sub-array split.
    boxes = [visualization.region_extent(region, frame=frame)]
    for ax, (label, pointed) in zip(axes, shown.items(), strict=True):
        visualization.probability_over_region(region, bands=bands, ax=ax,
                                              frame=frame, size=2.0)
        boxes.append(visualization.camera_rims(pointed, ax=ax, frame=frame,
                                               colors="type",
                                               type_names=("SST", "MST")))
        ax.set_title(label, loc="left")
        ax.set_xlabel("")
        ax.set_ylabel("")

    for ax in axes:
        visualization.frame_on(ax, *boxes)
    axes[-1].set_xlabel("along the region [deg]")
    axes[len(axes) // 2].set_ylabel("across it [deg]")

    # The panels keep an equal aspect, so a figure not shaped like their
    # contents shrinks the axes inside their slots rather than filling them.
    # The shape is only known once the widest configuration has been drawn, so
    # the figure is sized here rather than at construction.
    width = 9.6
    x_min, x_max, y_min, y_max = visualization._union(*boxes)
    panel = width * (y_max - y_min) / (x_max - x_min)
    figure.set_size_inches(width, len(axes) * (panel + 0.45) + 1.0)

    handles, labels = axes[0].get_legend_handles_labels()
    figure.tight_layout()
    legend = figure.legend(handles, labels, frameon=False, ncol=len(labels),
                           loc="lower center", bbox_to_anchor=(0.5, -0.01))
    for handle in legend.legend_handles:
        if hasattr(handle, "set_sizes"):
            handle.set_sizes([26])

    figure.savefig(args.figures / "footprints.png")
    plt.close(figure)


def coverage_figure(args, results):
    """Coverage against divergence, against tiles, and against sub-arrays."""
    figure, (left, right) = plt.subplots(1, 2, figsize=(9.6, 3.6))

    scan = results["strategies"]["scan"]
    left.plot([row["div"] for row in scan],
              [100 * row["covered"] for row in scan], marker="o", markersize=3,
              color="#2a78d6")
    best = results["strategies"]["divergent"]
    left.plot(best["div"], 100 * best["covered"], marker="*", markersize=14,
              color="#eb6834", linestyle="none")
    left.set_xlabel("divergence parameter $d$")
    left.set_ylabel("localization covered in stereo [%]")
    left.set_title("one divergent pointing", loc="left")
    left.set_ylim(0, 105)

    tiles = results["strategies"]["tiles"]
    right.plot(range(1, len(tiles) + 1), [100 * t["cumulative"] for t in tiles],
               marker="o", markersize=3, color="#1baf7a", label="sequential tiles")
    subarrays = results["strategies"]["subarray_scan"]
    right.plot([row["count"] for row in subarrays],
               [100 * row["covered"] for row in subarrays], marker="s",
               markersize=3, color="#eda100", label="simultaneous sub-arrays")
    right.axhline(100 * best["covered"], color="#eb6834", linestyle="--",
                  linewidth=1, label="best divergence")
    right.set_xlabel("pointings, or sub-arrays")
    right.set_title("moving the array, or splitting it", loc="left")
    right.set_ylim(0, 105)
    right.legend(frameon=False, loc="lower right")

    figure.savefig(args.figures / "coverage.png")
    plt.close(figure)


# -- writing it out --------------------------------------------------------

def substitutions(results):
    """Every number the prose quotes, as reStructuredText substitutions."""
    south, north = results["arrays"]["south"], results["arrays"]["north"]
    strat = results["strategies"]
    good = results["maps"]["bayestar_hlv"]
    first = results["maps"]["bayestar_hl"]
    prelim = results["maps"]["lalinference_prelim"]
    windows = results["visibility"]["sites"]

    south_window = windows["CTAO-South (Paranal)"]
    north_window = windows["CTAO-North (La Palma)"]

    values = {
        "gw-first-area": f"{first['area_deg2']:.0f} deg²",
        "gw-first-length": f"{first['max_separation_deg']:.0f}°",
        "gw-first-radius": f"{first['containment_radius_deg']:.1f}°",
        "gw-first-disc": f"{first['equivalent_radius_deg']:.1f}°",
        "gw-first-elongation": f"{first['elongation']:.0f}",
        "gw-good-area": f"{good['area_deg2']:.0f} deg²",
        "gw-good-radius": f"{good['containment_radius_deg']:.1f}°",
        "gw-good-delay": f"{good['since_merger_h']:.1f} h",
        "gw-prelim-radius": f"{prelim['containment_radius_deg']:.1f}°",
        "gw-radius-shrink":
            f"{first['containment_radius_deg'] / good['containment_radius_deg']:.0f}",
        "gw-area-shrink": f"{first['area_deg2'] / good['area_deg2']:.0f}",

        "gw-north-hours": f"{north_window['total_hours']:.0f} h",
        "gw-north-alt": f"{north_window['alt_max_deg']:.1f}°",
        "gw-south-hours": f"{south_window['total_hours']:.0f} h",
        "gw-south-alt": f"{south_window['alt_max_deg']:.0f}°",
        "gw-south-first": f"{south_window['first_h']:.1f} h",
        "gw-south-margin":
            f"{south_window['first_h'] - good['since_merger_h']:.1f} h",

        "gw-south-reach": f"{south['parallel_reach_deg']:.1f}°",
        "gw-north-reach": f"{north['parallel_reach_deg']:.1f}°",
        "gw-south-parallel-area": f"{south['parallel_area_deg2']:.0f} deg²",
        "gw-south-ceiling": f"{south['ceiling_div']:.3f}",
        "gw-south-ceiling-area": f"{south['ceiling_stereo_deg2']:.0f} deg²",
        "gw-south-half-omega": f"{south['half_omega_deg2']:.0f} deg²",
        "gw-south-ceiling-fraction": f"{100 * south['ceiling_fraction']:.0f} %",
        "gw-north-ceiling-fraction": f"{100 * north['ceiling_fraction']:.0f} %",
        "gw-subarray-advantage": f"{south['half_omega_deg2'] / south['ceiling_stereo_deg2']:.1f}",

        "gw-parallel-covered": f"{100 * strat['parallel']['covered']:.1f} %",
        "gw-div-covered": f"{100 * strat['divergent']['covered']:.1f} %",
        "gw-div-value": f"{strat['divergent']['div']:.3f}",
        "gw-div-spread": f"{strat['divergent']['spread_deg']:.0f}°",
        "gw-div-mean": f"{strat['divergent']['mean']:.1f}",
        "gw-div-on-target": str(strat["divergent"]["on_target"]),
        "gw-tiles": str(strat["tiles_to_match_divergence"]),
        "gw-tiles-99": str(strat["tiles_to_99"]),
        "gw-subarray-count": str(strat["subarrays"]["count"]),
        "gw-subarray-covered": f"{100 * strat['subarrays']['covered']:.1f} %",
        "gw-subarray-mean": f"{strat['subarrays']['mean']:.1f}",
        "gw-shaped-covered": f"{100 * strat['shaped']['covered']:.1f} %",
        "gw-shaped-mean": f"{strat['shaped']['mean']:.1f}",
        "gw-shaped-low": f"{strat['shaped']['profile'][0]:.1f}",
        "gw-shaped-high": f"{strat['shaped']['profile'][-1]:.1f}",
        "gw-shaped-mismatch": f"{strat['shaped']['mismatch']:.3f}",
        "gw-even-low": f"{strat['subarrays']['profile'][0]:.1f}",
        "gw-even-high": f"{strat['subarrays']['profile'][-1]:.1f}",
        "gw-weighted-count": str(strat["weighted"]["count"]),
        "gw-weighted-covered": f"{100 * strat['weighted']['covered']:.1f} %",
        "gw-weighted-high": f"{strat['weighted']['profile'][-1]:.1f}",
        "gw-weighted-mean": f"{strat['weighted']['mean']:.1f}",
        "gw-blur": f"{results['strategies']['camera_blur']:.3f}",
        "gw-telescopes": str(south["telescopes"]),
        "gw-region-area": f"{results['showcase_area_deg2']:.0f} deg²",
        "gw-region-dropped": f"{100 * results['showcase_dropped']:.0f} %",
        "gw-reach": f"{strat['shaped']['reach']:.0f}",
    }

    lines = [".. This file is generated by docs/scripts/make_gw170817.py.",
             ".. Every number the study page quotes lives here, so none of them",
             ".. can drift from the code that produced it. Do not edit.", ""]
    lines += [f".. |{name}| replace:: {value}" for name, value in values.items()]
    return "\n".join(lines) + "\n"


def _table(caption, header, rows, widths=None):
    """One `list-table` fragment."""
    out = [f".. list-table:: {caption}" if caption else ".. list-table::",
           "   :header-rows: 1"]
    if widths:
        out.append(f"   :widths: {' '.join(str(w) for w in widths)}")
    out.append("")
    for row in [header, *rows]:
        out.append(f"   * - {row[0]}")
        out += [f"     - {cell}" for cell in row[1:]]
    return "\n".join(out) + "\n"


def tables(results):
    """The page's tables, as fragments it includes."""
    maps_rows = []
    for name in ALERTS:
        entry = results["maps"][name]
        maps_rows.append([
            f"`{entry['label']} <https://gcn.nasa.gov/circulars/"
            f"{entry['gcn'].split()[-1]}>`_",
            f"+{entry['since_merger_h']:.1f} h",
            f"{entry['area_deg2']:.0f}",
            f"{entry['equivalent_radius_deg']:.1f}",
            f"**{entry['containment_radius_deg']:.1f}**",
            f"{entry['max_separation_deg']:.1f}",
            f"{entry['elongation']:.1f}",
        ])

    strat = results["strategies"]
    strategy_rows = [
        ["Parallel", f"{100 * strat['parallel']['covered']:.0f} %",
         f"{strat['parallel']['mean']:.0f}",
         f"{strat['parallel']['on_target']} of {strat['parallel']['telescopes']}",
         "1"],
        [f"Divergence, *d* = {strat['divergent']['div']:.3f}",
         f"{100 * strat['divergent']['covered']:.0f} %",
         f"{strat['divergent']['mean']:.0f}",
         f"{strat['divergent']['on_target']} of {strat['divergent']['telescopes']}",
         "1"],
        ["Sequential tiling", f"{100 * strat['tiling']['covered']:.0f} %",
         f"{strat['tiling']['mean']:.0f}",
         f"{strat['tiling']['on_target']} of {strat['tiling']['telescopes']}",
         str(strat["tiles_to_match_divergence"])],
        [f"Sub-arrays, {strat['subarrays']['count']} groups",
         f"{100 * strat['subarrays']['covered']:.0f} %",
         f"{strat['subarrays']['mean']:.0f}",
         f"{strat['subarrays']['on_target']} of {strat['subarrays']['telescopes']}",
         "1"],
        ["Shaped pointing", f"{100 * strat['shaped']['covered']:.0f} %",
         f"{strat['shaped']['mean']:.0f}",
         f"{strat['shaped']['on_target']} of {strat['shaped']['telescopes']}",
         "1"],
    ]

    def quintiles(row):
        return "  ".join(f"{value:.1f}" for value in row["profile"])

    # Two counts of each, so the reader can see that an even split stays flat
    # however many groups it is dealt into while a weighted one does not.
    unequal_rows = []
    for kind, scan, chosen in (("Even split", strat["subarray_scan"],
                                strat["subarrays"]["count"]),
                               ("Weighted split", strat["weighted_scan"],
                                strat["weighted"]["count"])):
        shown = sorted({chosen, max(row["count"] for row in scan)})
        for row in scan:
            if row["count"] in shown:
                unequal_rows.append([f"{kind}, {row['count']} groups",
                                     f"{100 * row['covered']:.1f} %",
                                     f"{row['mean']:.2f}", quintiles(row)])
    unequal_rows.append(["Shaped pointing",
                         f"{100 * strat['shaped']['covered']:.1f} %",
                         f"{strat['shaped']['mean']:.2f}",
                         quintiles(strat["shaped"])])

    beta_rows = []
    for beta in BETAS:
        row = strat["beta_scan"][f"{beta}"]
        mismatch = f"{row['mismatch']:.3f}"
        if row is min(strat["beta_scan"].values(),
                      key=lambda r: r["mismatch"]):
            mismatch = f"**{mismatch}**"
        beta_rows.append([f"{beta:.1f}", f"{100 * row['covered']:.1f} %",
                          f"{row['mean']:.1f}", mismatch,
                          f"{row['profile'][0]:.1f} → {row['profile'][-1]:.1f}"])

    parts = [
        ".. This file is generated by docs/scripts/make_gw170817.py. "
        "Do not edit.\n",
        ".. _gw170817-maps-table:\n",
        _table("The three sky maps circulated for GW170817, as measured here",
               ["Map", "Since merger", "90 % area [deg²]",
                "Equivalent radius [deg]", "Containment radius [deg]",
                "Length [deg]", "Elongation"],
               maps_rows, widths=[24, 10, 12, 14, 14, 10, 10]),
        "\n.. _gw170817-strategies-table:\n",
        _table("Five ways of covering the first GW170817 credible region from "
               "CTAO-South. Coverage is stereoscopic, and the mean is taken "
               "over the part covered",
               ["Strategy", "Coverage", "Telescopes per shower", "On target",
                "Exposures"],
               strategy_rows, widths=[30, 14, 18, 14, 12]),
        "\n.. _gw170817-unequal-table:\n",
        _table("Unequal sub-array splits against free pointing. Quintiles are "
               "mean multiplicity from the least to the most probable fifth",
               ["Strategy", "Coverage", "Mean", "Quintiles"],
               unequal_rows, widths=[30, 12, 10, 34]),
        "\n.. _gw170817-beta-table:\n",
        _table("Sweeping the tracking exponent β. Mismatch is the "
               "Kullback-Leibler divergence from the probability to the "
               "normalised multiplicity; zero would be exact proportionality",
               ["β", "Coverage", "Mean", "Mismatch", "Least to most likely fifth"],
               beta_rows, widths=[8, 12, 10, 12, 26]),
    ]
    return "\n".join(parts)


def parse_args(argv=None):
    here = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", type=Path, default=here / "_generated")
    parser.add_argument("--figures", type=Path,
                        default=here / "_static" / "gw170817")
    args = parser.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.figures.mkdir(parents=True, exist_ok=True)
    return args


def main(argv=None):
    args = parse_args(argv)

    region, dropped = placed(SHOWCASE)
    results = {
        "event": "GW170817",
        "merger_utc": MERGER,
        "site": SHOWCASE_SITE,
        "showcase": SHOWCASE,
        "credible_level": LEVEL,
        "m_cut": M_CUT,
        "placement": {"alt_deg": ALT.to_value(u.deg), "az_deg": AZ.to_value(u.deg)},
        "showcase_area_deg2": float(region.area.to_value(u.deg**2)),
        "showcase_dropped": dropped,
        "visibility": visibility(args),
        "maps": maps(args),
        "arrays": arrays(),
    }
    results["strategies"] = strategies(args, region,
                                       placed_bands(SHOWCASE))
    coverage_figure(args, results)

    (args.outdir / "gw170817_results.json").write_text(
        json.dumps(results, indent=2, default=str))
    (args.outdir / "gw170817_numbers.rst").write_text(substitutions(results))
    (args.outdir / "gw170817_tables.rst").write_text(tables(results))

    strat = results["strategies"]
    print(f"showcase: {len(region)} directions, "
          f"{results['showcase_area_deg2']:.0f} deg2 above the horizon")
    for name in ("parallel", "divergent", "tiling", "subarrays", "shaped"):
        row = strat[name]
        print(f"  {name:<11} {100 * row['covered']:5.1f}%  "
              f"<m> {row['mean']:5.2f}  {row['on_target']:>2}/{row['telescopes']} "
              f"on target, {row['exposures']} exposure(s)")


if __name__ == "__main__":
    main()

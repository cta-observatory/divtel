#!/usr/bin/env python
"""
Precompute the GW170817 credible regions that ship with divtel.

Reading a published HEALPix localization needs ``astropy-healpix`` (the
``divtel[skymap]`` extra) and a download from the LIGO document server. Neither
is available where these regions are most wanted: the documentation build, and
the marimo notebooks that run in a reader's browser under Pyodide, where the
extension module has no wheel and the document server will not answer a
cross-origin request.

So the regions are cut once, here, and committed as plain tables of right
ascension, declination and probability. `divtel.region.SkyRegion.from_table`
reads them back, and from there every strategy in `divtel.strategy` works with
no sky-map machinery at all.

Run it from a checkout when the maps or the credible levels need to change::

    pip install -e .[skymap]
    python divtel/data/gw170817/make_regions.py

Sky maps are downloaded into a cache directory (``--cache``) and kept, so a
rerun costs nothing. ``--no-download`` insists on the cache and fails rather
than reaching for the network.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from divtel import skymap

HERE = Path(__file__).resolve().parent

# The alert sequence as observers received it. The final catalogue map is left
# out on purpose: a follow-up is planned against the best localization available
# when the telescope could move, and nobody had GWTC-1 in 2017.
MAPS = ["bayestar_hl", "bayestar_hlv", "lalinference_prelim"]

# 0.9 is the region everyone quotes and the one the study is scored on. 0.5 is
# the core, for drawing. 0.95 costs little and lets a reader ask what one more
# credible point buys, which on an arc is a great deal.
LEVELS = [0.5, 0.9, 0.95]

# nside 256 is 0.052 deg**2 per pixel: two orders of magnitude finer than the
# cameras that have to cover these regions, and small enough that the largest
# region here is a few thousand rows.
NSIDE = 256


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cache", type=Path, default=HERE / "_maps",
                        help="where to keep the downloaded FITS maps")
    parser.add_argument("--outdir", type=Path, default=HERE)
    parser.add_argument("--nside", type=int, default=NSIDE)
    parser.add_argument("--no-download", action="store_true",
                        help="use the cache only, and fail if a map is missing")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    for name in MAPS:
        path = skymap.fetch(name, args.cache, download=not args.no_download)
        sky_map = skymap.load(path, skymap.SKYMAPS[name], max_nside=args.nside)

        for level in LEVELS:
            region = skymap.region(sky_map, level)
            region.meta["nside"] = int(sky_map.nside)
            region.meta["source"] = skymap.SKYMAPS[name].url

            out = args.outdir / f"{name}_{round(100 * level)}.ecsv.gz"
            region.to_table(out)
            print(f"{out.name:<34} {len(region):>6} pixels  "
                  f"{region.area.value:7.1f} deg2  "
                  f"{region.weights.sum():.4f} probability")


if __name__ == "__main__":
    main()

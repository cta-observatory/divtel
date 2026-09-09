# Changelog

## Unreleased

**Pointing a whole region of sky**

divtel could point an array at a direction. It can now cover a *region* with
structure — a survey field, a galaxy catalogue, a gravitational-wave
localization — by three strategies that are not divergence.

- `divtel.region.SkyRegion`: a set of sky directions with a weight on each,
  which every strategy works on. Places itself over an array at a chosen
  elevation, cuts itself at the horizon, and measures its own size three ways
  (equivalent radius, containment radius, length — for an arc those differ by
  an order of magnitude, and only the second maps onto what an array must
  span).
- `divtel.strategy`: sequential tiling (`tile_region`), simultaneous
  sub-arrays (`point_subarrays`, `point_subarrays_by_type`, and
  `weighted_split` for groups sized by the probability they catch), and
  shaped pointing (`shaped_pointing`, `shaped_pointing_by_type`) that sets
  every telescope independently so the array's depth follows the weight.
  Plus `reach`, `describe`, `camera_blur` and `multiplicity_by_probability`
  to score any of them.
- `divtel.pointing` gains the divergence solvers: `pointing_spread`,
  `div_for_half_angle`, `div_for_multiplicity`, `best_pointing`, `div_scan`,
  `spread_scan` and `solve_div`.
- `divtel.observation` gains `observable_windows`, `altaz_track` and
  `Window`: when a target is up and the sky is dark.
- `Array.group_by` accepts `"fov_radius"` beside `"camera_radius"`. Angular
  radius is what the strategies mean by a telescope's type.
- `divtel.visualization` gains `camera_rims` (every telescope's field of
  view drawn over the region, sampled on the sphere so a camera at the edge
  of the frame keeps its true shape), `probability_over_region`,
  `multiplicity_over_region`, `multiplicity_by_probability`, `sky_bands`,
  and the `Projection` frame they share — `projection_frame`,
  `project_region`, `project_directions`, `region_extent`, `region_span`,
  `frame_on` and `type_colors`.

**Sky maps** (optional: `pip install divtel[skymap]`)
- `divtel.skymap` reads a published HEALPix localization, cuts credible
  regions out of it, splits a region into its disconnected lobes, and hands
  the result on as a `SkyRegion`. Needs `astropy-healpix`, which nothing else
  in divtel does.

**Data**
- The three GW170817 alert maps ship as precomputed credible regions in
  `divtel/data/gw170817`, at three credible levels each. They need neither
  the extra nor the network, which is what lets the browser notebooks use
  them. `make_regions.py` beside them is the script that cut them.

**Docs**
- New study: [Covering GW170817](https://cta-observatory.github.io/divtel/gw170817.html).
  Divergence reaches 70 % of the first alert map and can do no better at any
  value; sub-arrays cover all of it in one exposure; shaped pointing then
  makes multiplicity track the probability. Every number on the page is
  computed at build time by `docs/scripts/make_gw170817.py`, so none of them
  can drift from the code.
- New interactive notebook, `gw170817_strategies`, embedded in that page.

## v1.1.0 (2026-09-08)

**Hyper field of view**
- `hyper_fov`, `display_hyper_fov` and `multiplicity_plot` rename their
  `m_cut` argument to `min_telescopes` and default it to 2 instead of 1: sky
  seen by only one telescope isn't stereo-reconstructable, so it no longer
  counts by default. Pass `min_telescopes=1` to get the old behavior.

**Data**
- Fixed the South array (`cta-south-paranal-alpha-prod6.ecsv`) telescope
  specs, which were wrong.

**Docs**
- Sub-arrays marimo notebook now covers both CTAO sites.
- New 3D divergence-geometry marimo notebook, wired into
  `interactive_display`.
- `CITATION.cff` added, derived from `codemeta.json`.

## v1.0.0 (2026-08-28)

First stable release. Everything below shipped since `v0.1` (April 2022).

**Layouts**
- Load arrays from ECSV files (`divtel.layout.load_array`,
  `divtel.layout.load_table`), including a camera radius given as either a
  length or an angle on the sky.
- CTAO's `id`-based telescope numbering is now semantic and used by
  `Array.group_by`.

**Pointing**
- `div` is validated to `[0, 1]` and documented precisely (a fixed 100 m
  reference distance, not a raw angle).
- Divergent pointing clamps at the horizon instead of pointing telescopes
  through the ground on sloped arrays.
- A single telescope can be pointed at an object or directly in alt/az.
- `Telescope.position` and the rest of the geometry API are now unit-aware
  throughout (`astropy.units.Quantity`).

**Hyper field of view**
- `Array.hyper_fov` reworked onto a proper equal-area sky projection, fixing
  a real discontinuity near the map boundary.
- `Array.multiplicity_profile` and `Array.multiplicity_moments` report how
  *well*, not just how *much*, sky is covered.

**Sub-arrays**
- `Array.group_by` splits a mixed array (e.g. LSTs vs MSTs) into independent
  sub-arrays that share state with the parent.

**Observation**
- New `divtel.observation.Observation` ties ground-frame pointing to a real
  site and time, converting sky coordinates to alt/az and back.

**Export**
- `Array.export_cfg` writes a `sim_telarray` configuration file for a
  pointed array.

**Docs**
- Full user guide and API reference, plus an in-browser interactive demo
  (a marimo notebook compiled to WebAssembly) published via GitHub Pages.
- Three worked tutorial notebooks.

**Packaging**
- Migrated to PEP 621 `pyproject.toml` with `setuptools_scm` for
  git-derived versioning, PEP 735 dependency groups, and `uv` for local
  development.

**Fixes**
- `random_array` notebook helper no longer drifts off-center (#14).
- Azimuth sign in `Telescope.point_to_object`.
- Various matplotlib compatibility and CI hygiene fixes.

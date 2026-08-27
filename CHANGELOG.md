# Changelog

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

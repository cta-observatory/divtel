# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import json
import os
import subprocess
import sys
import shutil
from pathlib import Path

sys.path.insert(0, os.path.abspath('..'))
import divtel

HERE = Path(__file__).parent.resolve()
REPO_ROOT = HERE.parent

# -- Project information -----------------------------------------------------

project = 'divtel'
copyright = '2022, Thomas Vuillaume'
author = 'Thomas Vuillaume'

# The full version, including alpha/beta/rc tags
release = divtel.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinxarg.ext',
    'sphinx.ext.napoleon',
    'myst_parser',
]

# A dataclass's fields are documented twice otherwise: once by autodoc from the
# annotations, once by napoleon from the Attributes section. Rendering the
# section as :ivar: fields on the class leaves one description of each.
napoleon_use_ivar = True


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '_generated']

source_suffix = [
    '.rst',
    '.md',
]

# Figures get numbered captions, and :numref: in the prose resolves to
# "Figure N" rather than a bare "the figure below".
numfig = True

# The generated numbers are substitutions, and a substitution has to be defined
# in the document that uses it. Putting them in the prolog defines them once for
# every page rather than making each page include them.
rst_prolog = '.. include:: /_generated/gw170817_numbers.rst\n'

# -- Interactive (marimo / WebAssembly) --------------------------------------
# GitHub Pages only serves static files, so the interactive demo is a marimo
# notebook exported to WebAssembly: it carries its own Python and runs entirely
# in the reader's browser.

MARIMO_DIR = REPO_ROOT / 'examples' / 'marimo'

# Exported in this order; each becomes marimo/<name>/index.html on the site.
MARIMO_NOTEBOOKS = [
    'interactive_display',
    'two_sites',
    'observing_a_source',
    'choosing_div',
    'gw170817_strategies',
]

# -- Static study plots -------------------------------------------------------
# The ceiling and studies pages illustrate the argument with matplotlib
# figures that are not the interactive marimo demos above. They are rendered
# here rather than committed, for the same reason the gw170817 figures are:
# a plot in the repository can silently drift from the code that made it.

STATIC_PLOTS = [
    HERE / '_static' / 'studies' / 'array_and_sky.png',
    HERE / '_static' / 'studies' / 'div_schema.png',
    HERE / '_static' / 'studies' / 'tracking.png',
    HERE / '_static' / 'ceiling' / 'fig3_north.png',
    HERE / '_static' / 'ceiling' / 'fig3_south.png',
    HERE / '_static' / 'ceiling' / 'ceiling_maps.png',
]


def _sweep(array, groups, divs, alts):
    """Coverage and mean multiplicity against divergence, whole array and by
    telescope type, at each of ``alts``."""
    import numpy as np
    import astropy.units as u

    whole, by_type = {}, {}
    for alt in alts:
        areas, means = [], []
        for div in divs:
            array.divergent_pointing(div, alt, 180 * u.deg)
            area, patches = array.hyper_fov(min_telescopes=2)
            areas.append(area.to_value(u.deg**2))
            means.append(array.multiplicity_moments(patches=patches)[0])
        whole[int(alt.to_value(u.deg))] = (np.array(areas), np.array(means))

    for name, subset in array.group_by(groups).items():
        areas, means = [], []
        for div in divs:
            subset.divergent_pointing(div, 90 * u.deg, 180 * u.deg)
            area, patches = subset.hyper_fov(min_telescopes=2)
            areas.append(area.to_value(u.deg**2))
            means.append(subset.multiplicity_moments(patches=patches)[0])
        by_type[name] = (np.array(areas), np.array(means))
    return whole, by_type


def _find_ceiling_div(divs, means):
    """The first divergence at which mean multiplicity drops to 2."""
    import numpy as np

    below = np.where(means <= 2)[0]
    if len(below) == 0:
        return divs[-1]
    i = below[0]
    if i == 0:
        return divs[0]
    x0, x1 = divs[i - 1], divs[i]
    y0, y1 = means[i - 1], means[i]
    return float(x0 + (2 - y0) * (x1 - x0) / (y1 - y0))


def _ceiling_div(array):
    """The array's own ceiling divergence, at 60 degrees elevation."""
    import numpy as np
    import astropy.units as u

    divs = np.linspace(0.001, 0.15, 80)
    means = []
    for div in divs:
        array.divergent_pointing(div, 60 * u.deg, 180 * u.deg)
        _, patches = array.hyper_fov(min_telescopes=2)
        means.append(array.multiplicity_moments(patches=patches)[0])
    return _find_ceiling_div(divs, np.array(means))


def _generate_static_plots(app):
    """Render the ceiling and studies pages' figures before Sphinx reads them."""
    if app.builder.name != 'html':
        return

    import numpy as np
    import astropy.units as u
    import matplotlib.pyplot as plt
    from importlib.resources import files
    from astropy.coordinates import AltAz, SkyCoord, get_body
    from astropy.utils import iers

    from divtel.layout import load_array
    from divtel.observation import Observation
    from divtel.visualization import display_hyper_fov
    from divtel.telescope import Telescope, Array
    from divtel.pointing import pointG_position

    for path in STATIC_PLOTS:
        path.parent.mkdir(parents=True, exist_ok=True)

    array_n = load_array(files('divtel') / 'data' / 'cta-north-lapalma-alpha-prod6.ecsv')
    array_s = load_array(files('divtel') / 'data' / 'cta-south-paranal-alpha-prod6.ecsv')

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    array_n.divergent_pointing(0.03, 70 * u.deg, 180 * u.deg)
    array_n.display_2d(projection='xy', ax=axes[0])
    axes[0].set_title('CTAO-North on the ground')
    axes[0].legend(loc='lower right')
    display_hyper_fov(array_n, ax=axes[1])
    fig.tight_layout()
    fig.savefig(HERE / '_static' / 'studies' / 'array_and_sky.png', dpi=200)
    plt.close(fig)

    # A toy array whose divergent pointing is computed the same way as the
    # real arrays above, so the schematic is a faithful picture of the
    # geometry `pointG_position` / `tel_div_pointing` actually produce,
    # rather than a hand-drawn stand-in.
    div = 0.12
    tel_x = np.array([-300.0, -180.0, -60.0, 60.0, 180.0, 300.0]) * u.m
    toy = Array([Telescope(x, 0 * u.m, 0 * u.m, focal=16 * u.m,
                            camera_radius=1 * u.m, tel_id=i)
                 for i, x in enumerate(tel_x, start=1)])
    toy.divergent_pointing(div, 90 * u.deg, 180 * u.deg)

    positions = toy.positions_array.to_value(u.m)
    vectors = toy.pointing_vectors
    barycentre = toy.barycenter.to_value(u.m)
    g_point = pointG_position(toy.barycenter, div, 90 * u.deg,
                               180 * u.deg).to_value(u.m)

    fig, ax = plt.subplots(figsize=(5.5, 8))
    arrow_len = 220
    for (x, _, z), (vx, _, vz) in zip(positions, vectors):
        ax.plot([g_point[0], x], [g_point[2], z], '--', color='tab:orange', lw=1)
        ax.arrow(x, z, vx * arrow_len, vz * arrow_len,
                 width=3, head_width=16, head_length=22, color='k',
                 length_includes_head=True, zorder=3)

    ax.scatter(positions[:, 0], positions[:, 2], color='tab:blue', s=40, zorder=4)
    ax.annotate('telescope', xy=(positions[-1, 0], positions[-1, 2]),
                xytext=(12, 8), textcoords='offset points')
    ax.scatter([barycentre[0]], [barycentre[2]], color='black', s=25, zorder=4)
    ax.annotate('B', xy=(barycentre[0], barycentre[2]),
                xytext=(barycentre[0] + 12, barycentre[2] + 10))
    ax.scatter([g_point[0]], [g_point[2]], color='tab:orange', marker='*',
               s=160, zorder=4)
    ax.annotate('G', xy=(g_point[0], g_point[2]),
                xytext=(g_point[0] + 15, g_point[2] - 5))

    # |BG|: the barycenter-to-G distance the divergence angle is measured against.
    ax.annotate('', xy=(barycentre[0], g_point[2]), xytext=(barycentre[0], barycentre[2]),
                arrowprops=dict(arrowstyle='-', lw=1.3, color='gray'))
    ax.annotate('|BG|', xy=(barycentre[0] - 15, (barycentre[2] + g_point[2]) / 2),
                ha='right', va='center', color='gray')

    # D and theta_D: one telescope's baseline from the barycenter and the
    # angle its pointing diverges by, i.e. div = sin(theta_D) at D = 100 m.
    i_ref = 4
    x_ref, z_ref = positions[i_ref, 0], positions[i_ref, 2]
    vx_ref, vz_ref = vectors[i_ref, 0], vectors[i_ref, 2]
    ax.annotate('', xy=(x_ref, -25), xytext=(barycentre[0], -25),
                arrowprops=dict(arrowstyle='<->', lw=1))
    ax.annotate('D', xy=((barycentre[0] + x_ref) / 2, -55), ha='center')

    ax.plot([x_ref, x_ref], [z_ref, z_ref + 200], color='gray', linestyle=':', lw=1.2)
    arc_r = 60
    angles = np.linspace(np.pi / 2, np.arctan2(vz_ref, vx_ref), 20)
    ax.plot(x_ref + arc_r * np.cos(angles), z_ref + arc_r * np.sin(angles),
            color='k', lw=1)
    ax.annotate(r'$\theta_D$', xy=(x_ref + arc_r * 1.5, z_ref + arc_r * 1.1))

    ax.set_xlim(-420, 420)
    ax.set_ylim(-1050, 300)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('ground position [m]')
    ax.set_ylabel('along the mean pointing [m]')
    ax.set_title(f'Divergent pointing, div = {div:.2f}')
    fig.tight_layout()
    fig.savefig(HERE / '_static' / 'studies' / 'div_schema.png', dpi=200)
    plt.close(fig)

    divs = np.linspace(0.002, 0.16, 40)
    alts = np.array([30, 50, 70, 90]) * u.deg
    whole_n, by_type_n = _sweep(array_n, {'LST': range(1, 5), 'MST': range(5, 14)},
                                divs, alts)
    whole_s, by_type_s = _sweep(array_s, {'MST': range(1, 15), 'SST': range(15, 52)},
                                divs, alts)

    for whole, by_type, site, path in (
        (whole_n, by_type_n, 'CTAO-North', 'fig3_north.png'),
        (whole_s, by_type_s, 'CTAO-South', 'fig3_south.png'),
    ):
        fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.5))
        for alt, (areas, means) in whole.items():
            line, = left.plot(divs, areas, label=f'{alt}° (whole array)')
            right.plot(divs, means, color=line.get_color())
        for name, (areas, means) in by_type.items():
            line, = left.plot(divs, areas, '--', label=f'{name} (zenith)')
            right.plot(divs, means, '--', color=line.get_color())
        right.axhline(2, color='k', linewidth=1, linestyle=':',
                      label='stereoscopic floor')
        left.set_xlabel('div')
        left.set_ylabel('stereoscopic hyper FoV [deg$^2$]')
        right.set_xlabel('div')
        right.set_ylabel('mean multiplicity')
        left.grid(alpha=0.3)
        right.grid(alpha=0.3)
        left.legend(frameon=False, fontsize=8)
        right.legend(frameon=False, fontsize=8)
        fig.suptitle(f'{site}: stereoscopic coverage and mean multiplicity '
                     'against div')
        fig.tight_layout()
        fig.savefig(HERE / '_static' / 'ceiling' / path, dpi=200)
        plt.close(fig)

    div_n = _ceiling_div(array_n)
    div_s = _ceiling_div(array_s)
    array_n.divergent_pointing(div_n, 60 * u.deg, 180 * u.deg)
    array_s.divergent_pointing(div_s, 60 * u.deg, 180 * u.deg)
    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 5))
    display_hyper_fov(array_n, ax=left)
    left.set_title(f'CTAO-North, div = {div_n:.3f}')
    display_hyper_fov(array_s, ax=right)
    right.set_title(f'CTAO-South, div = {div_s:.3f}')
    fig.suptitle('Coverage at the ceiling divergence, shaded by multiplicity')
    fig.tight_layout()
    fig.savefig(HERE / '_static' / 'ceiling' / 'ceiling_maps.png', dpi=200)
    plt.close(fig)

    # IERS tables would otherwise be fetched over the network at build time,
    # which CI cannot rely on; the bundled tables are close enough for a
    # single illustrative night.
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    target = SkyCoord(ra=83.633 * u.deg, dec=22.015 * u.deg)
    start = Observation(site='north', time='2026-12-01T00:00:00')
    coarse = start.time + np.linspace(-12, 12, 49) * u.hour
    sun = get_body('sun', coarse, location=start.location).transform_to(
        AltAz(obstime=coarse, location=start.location))
    midnight = Observation(site='north', time=coarse[sun.alt.argmin()])
    hours = np.arange(-6, 6.01, 0.5)
    times, areas, means = [], [], []
    for offset in hours:
        moment = midnight.after(offset * u.hour)
        alt, az = moment.altaz_of(target)
        if alt <= 0 * u.deg:
            continue
        array_n.divergent_pointing(0.04, alt, az)
        times.append(offset)
        areas.append(array_n.hyper_fov()[0].to_value(u.deg**2))
        means.append(array_n.multiplicity_moments()[0])
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(times, areas, marker='o', color='tab:blue')
    ax.set_xlabel('hours from the middle of the night')
    ax.set_ylabel('hyper FoV [deg$^2$]', color='tab:blue')
    ax.tick_params(axis='y', labelcolor='tab:blue')
    twin = ax.twinx()
    twin.plot(times, means, marker='s', color='tab:red')
    twin.set_ylabel('mean multiplicity', color='tab:red')
    twin.tick_params(axis='y', labelcolor='tab:red')
    ax.set_title('Tracking Crab Nebula at div = 0.04')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / '_static' / 'studies' / 'tracking.png', dpi=200)
    plt.close(fig)


def _cleanup_static_plots(app, exception):
    """Drop the rendered figures from the source tree once the build is done."""
    for path in STATIC_PLOTS:
        path.unlink(missing_ok=True)


def _build_divtel_wheel(destination):
    """Build a divtel wheel for the in-browser interpreter.

    Pyodide installs from wheels, and the divtel release on PyPI is an
    outdated sdist, so the wheel is built from the working tree and published
    alongside the exported notebook.
    """
    # uv-managed virtualenvs have no pip, and pip-managed ones have no uv, so
    # try both rather than pinning the docs build to one workflow.
    builders = []
    if shutil.which('uv'):
        builders.append(['uv', 'build', '--wheel', '--out-dir', str(destination),
                         str(REPO_ROOT)])
    builders.append([sys.executable, '-m', 'pip', 'wheel', '--no-deps',
                     '--wheel-dir', str(destination), str(REPO_ROOT)])

    for cmd in builders:
        if subprocess.run(cmd).returncode == 0:
            break
    else:
        raise RuntimeError('could not build a divtel wheel for the marimo export')

    # `uv build` drops a .gitignore of "*" beside its output, on the assumption
    # that build artefacts are never committed. Here they are the payload, and
    # publishing to gh-pages goes through git: left in place it hides the wheel
    # from the deploy, and the published notebook cannot install divtel.
    (destination / '.gitignore').unlink(missing_ok=True)

    return sorted(destination.glob('divtel-*.whl'))[0]


def _export_one_marimo(notebook, outdir, wheel):
    """Export a single marimo notebook to WebAssembly, ready to serve."""
    if outdir.exists():
        shutil.rmtree(outdir)

    export = subprocess.run(
        [sys.executable, '-m', 'marimo', 'export', 'html-wasm',
         str(notebook), '-o', str(outdir), '--mode', 'run'],
        capture_output=True,
        text=True,
    )
    if export.returncode != 0:
        # Sphinx reports only the exit status of a failed handler, which says
        # nothing about what marimo objected to. Carry its output into the
        # exception so the build log explains itself.
        raise RuntimeError(
            f'marimo export of {notebook.name} failed (exit {export.returncode}):\n'
            f'{export.stdout}\n{export.stderr}'
        )

    # marimo seeds its exports with editor scaffolding; it has no business
    # being published as part of the documentation.
    (outdir / 'CLAUDE.md').unlink(missing_ok=True)

    # Each notebook resolves the wheel relative to its own page, so every
    # export carries its own copy. The wheel is a few tens of kilobytes, which
    # is far cheaper than teaching the notebooks how deep they are published.
    #
    # Do not rename this to `wheels`: publishing to gh-pages honours the
    # repository's .gitignore, and the standard Python template ignores
    # `wheels/`. The directory built fine and was then dropped on the way out,
    # leaving the published notebook unable to install divtel.
    wheelhouse = outdir / 'pypi'
    wheelhouse.mkdir(parents=True, exist_ok=True)
    shutil.copy2(wheel, wheelhouse / wheel.name)
    # setuptools_scm bakes the version into the wheel filename, so record it in
    # a manifest the notebook can read instead of hardcoding it there.
    (wheelhouse / 'manifest.json').write_text(json.dumps({'divtel': wheel.name}))


def _export_marimo(app, exception):
    """Export every marimo notebook to WebAssembly into the built site."""
    if exception is not None or app.builder.name != 'html':
        return

    # marimo shells out to uv to resolve the notebook's imports, so a machine
    # without it fails here rather than at install time.
    if not shutil.which('uv'):
        raise RuntimeError(
            'uv is required to export the marimo notebooks to WebAssembly; '
            'install it with `pip install uv`'
        )

    root = Path(app.outdir) / 'marimo'
    if root.exists():
        shutil.rmtree(root)

    # Built once and copied into each export, rather than once per notebook:
    # the wheel is identical every time and building it is the slow part.
    staging = root / '_wheel'
    staging.mkdir(parents=True, exist_ok=True)
    wheel = _build_divtel_wheel(staging)

    for name in MARIMO_NOTEBOOKS:
        notebook = MARIMO_DIR / f'{name}.py'
        if not notebook.exists():
            raise RuntimeError(f'{notebook} is listed in MARIMO_NOTEBOOKS but missing')
        _export_one_marimo(notebook, root / name, wheel)

    shutil.rmtree(staging)


# -- Generated study material -------------------------------------------------
# The GW170817 page quotes a great many numbers and must not be able to disagree
# with the code. None of them is typed into the page: this script computes them,
# writes the figures, and writes the tables and substitutions the page includes.
# It runs before Sphinx reads a source file, so the fragments exist by the time
# the page is parsed.

GW170817_SCRIPT = HERE / 'scripts' / 'make_gw170817.py'


def _make_gw170817(app):
    """Compute the GW170817 study's figures, tables and numbers."""
    if os.environ.get('DIVTEL_DOCS_SKIP_GW170817'):
        # Worth setting while editing prose. Not worth setting before
        # publishing: the page then carries whatever the last run left behind.
        return

    result = subprocess.run([sys.executable, str(GW170817_SCRIPT)],
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f'{GW170817_SCRIPT.name} failed (exit {result.returncode}):\n'
            f'{result.stdout}\n{result.stderr}'
        )
    print(result.stdout, end='')


def setup(app):
    app.connect('builder-inited', _make_gw170817)
    app.connect('builder-inited', _generate_static_plots)
    app.connect('build-finished', _cleanup_static_plots)
    app.connect('build-finished', _export_marimo)


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'furo'

html_theme_options = {
    'sidebar_hide_name': True,
}

html_static_path = ['_static']
html_css_files = ['custom.css']

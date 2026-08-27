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


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

source_suffix = [
    '.rst',
    '.md',
]

# -- Interactive (marimo / WebAssembly) --------------------------------------
# GitHub Pages only serves static files, so the interactive demo is a marimo
# notebook exported to WebAssembly: it carries its own Python and runs entirely
# in the reader's browser.

MARIMO_DIR = REPO_ROOT / 'examples' / 'marimo'

# Exported in this order; each becomes marimo/<name>/index.html on the site.
MARIMO_NOTEBOOKS = [
    'interactive_display',
    'sub_arrays',
    'observing_a_source',
    'choosing_div',
]


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


def setup(app):
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

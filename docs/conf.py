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

import divtel
sys.path.insert(0, os.path.abspath('..'))

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

MARIMO_NOTEBOOK = REPO_ROOT / 'examples' / 'marimo' / 'interactive_display.py'


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

    return sorted(destination.glob('divtel-*.whl'))[0]


def _export_marimo(app, exception):
    """Export the marimo notebook to WebAssembly into the built site."""
    if exception is not None or app.builder.name != 'html':
        return

    outdir = Path(app.outdir) / 'marimo'
    if outdir.exists():
        shutil.rmtree(outdir)

    # marimo shells out to uv to resolve the notebook's imports, so a machine
    # without it fails here rather than at install time.
    if not shutil.which('uv'):
        raise RuntimeError(
            'uv is required to export the marimo notebook to WebAssembly; '
            'install it with `pip install uv`'
        )

    export = subprocess.run(
        [sys.executable, '-m', 'marimo', 'export', 'html-wasm',
         str(MARIMO_NOTEBOOK), '-o', str(outdir), '--mode', 'run'],
        capture_output=True,
        text=True,
    )
    if export.returncode != 0:
        # Sphinx reports only the exit status of a failed handler, which says
        # nothing about what marimo objected to. Carry its output into the
        # exception so the build log explains itself.
        raise RuntimeError(
            f'marimo export failed (exit {export.returncode}):\n'
            f'{export.stdout}\n{export.stderr}'
        )

    # marimo seeds its exports with editor scaffolding; it has no business
    # being published as part of the documentation.
    (outdir / 'CLAUDE.md').unlink(missing_ok=True)

    # setuptools_scm bakes the version into the wheel filename, so record it in
    # a manifest the notebook can read instead of hardcoding it there.
    wheelhouse = outdir / 'wheels'
    wheelhouse.mkdir(parents=True, exist_ok=True)
    wheel = _build_divtel_wheel(wheelhouse)
    (wheelhouse / 'manifest.json').write_text(json.dumps({'divtel': wheel.name}))


def setup(app):
    app.connect('build-finished', _export_marimo)


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# No custom static files yet; point html_static_path at a directory when there
# are some, otherwise Sphinx warns about the missing path on every build.
html_static_path = []

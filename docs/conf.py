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

# Notebooks
notebook_dir = '../examples/notebooks/'
os.makedirs('examples', exist_ok=True)
[shutil.copy(notebook_dir + file, 'examples') for file in os.listdir(notebook_dir) if file.endswith('.ipynb')]


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
    'nbsphinx',
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

# nbsphinx
nbsphinx_allow_errors = False
nbsphinx_execute = 'always'  # disable with 'never', force with 'always'

# Point readers at the live version. The rendered notebook shows ipywidgets
# sliders, but they cannot move: `@interact` re-runs Python on every change and
# a static host has no kernel to run it. The marimo build below does.
nbsphinx_prolog = """
.. raw:: html

    <div class="admonition tip">
      <p class="admonition-title">Interactive version</p>
      <p>
        The sliders below are a static snapshot &mdash; moving them needs a
        Python kernel, which a static page has not got.
        <a href="../marimo/index.html"><strong>Open the interactive
        version</strong></a> for sliders that work. Nothing is installed and
        nothing is uploaded: Python runs on your own machine, compiled to
        WebAssembly.
      </p>
    </div>
"""


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

    subprocess.run(
        [sys.executable, '-m', 'marimo', 'export', 'html-wasm',
         str(MARIMO_NOTEBOOK), '-o', str(outdir), '--mode', 'run'],
        check=True,
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

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

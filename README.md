# divtel: Divergent pointing mode for Imaging Atmospheric Cherenkov Telescopes arrays


<div align="left">

[![Build status](https://github.com/cta-observatory/divtel/workflows/build/badge.svg?branch=master&event=push)](https://github.com/cta-observatory/divtel/actions?query=workflow%3Abuild)
[![Python Version](https://img.shields.io/pypi/pyversions/divtel.svg)](https://pypi.org/project/divtel/)
[![Semantic Versions](https://img.shields.io/badge/%20%20%F0%9F%93%A6%F0%9F%9A%80-semantic--versions-e10079.svg)](https://github.com/cta-observatory/divtel/releases)
[![License](https://img.shields.io/github/license/cta-observatory/divtel?style=flat)](https://github.com/cta-observatory/divtel/blob/master/LICENSE)
[![Documentation](https://img.shields.io/github/workflow/status/cta-observatory/divtel/Sphinx%20docs%20to%20gh-pages/master?label=Documentation)](https://cta-observatory.github.io/divtel/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6415137.svg)](https://doi.org/10.5281/zenodo.6415137)
</div>

divtel makes toy simulations for the divergent pointing mode of Imaging Atmospheric Cherenkov Telescopes arrays.

## 👨‍💻 Install

### Users
``` 
pip install divtel
```

### Developers

With [uv](https://docs.astral.sh/uv/) (recommended):

```
git clone https://github.com/cta-observatory/divtel.git
cd divtel
uv sync
```

`uv sync` creates a virtual environment in `.venv`, installs divtel in editable
mode, and pulls in the development dependencies (pytest, sphinx, ruff) declared
as [PEP 735](https://peps.python.org/pep-0735/) dependency groups. Add
`--extra examples` if you also want to run the notebooks in `examples/`.

With pip (requires pip >= 25.1 for `--group`):

```
git clone https://github.com/cta-observatory/divtel.git
cd divtel
pip install -e . --group dev
```

Then run the tests:

```
pytest
```

> **Note:** install divtel before importing it, even from a source checkout.
> The version is derived from git by
> [setuptools_scm](https://setuptools-scm.readthedocs.io/) at install time and
> written to `divtel/_version.py`; importing an uninstalled source tree reports
> `__version__ == "0.0.0"`.

## 📚 Documentation

The docs are published to <https://cta-observatory.github.io/divtel/>. To build
them locally:

```
uv sync --group docs --extra examples
sphinx-build -b html docs docs/_build/html
```

GitHub Pages serves static files only, so it cannot preview the result:
`file://` will not work, and the interactive demo needs a real HTTP origin.
Serve the build instead:

```
python -m http.server 8000 -d docs/_build/html
```

### The interactive demo

`examples/marimo/interactive_display.py` is a [marimo](https://marimo.io)
notebook. Sphinx exports it to WebAssembly during the build, so the published
page ships its own Python interpreter and runs entirely in the reader's
browser — sliders included, with no server and nothing to install.

That export shells out to [uv](https://docs.astral.sh/uv/), so `uv` must be on
your `PATH` to build the docs — `uv sync` above already provides it, but a
pip-based checkout needs `pip install uv`.


To work on the demo:

```
marimo edit examples/marimo/interactive_display.py
```


## 🛡 License

[![License](https://img.shields.io/github/license/cta-observatory/divtel?style=flat)](https://github.com/cta-observatory/divtel/blob/master/LICENSE)

This project is licensed under the terms of the `MIT` license. See [LICENSE](https://github.com/cta-observatory/divtel/blob/master/LICENSE) for more details.

## 📃 Citation

```bibtex
@software{thomas_vuillaume_2022_6415138,
  author       = {Thomas Vuillaume and
                  Alice Donini and
                  Thomas Gasparetto},
  title        = {cta-observatory/divtel: v0.1},
  month        = apr,
  year         = 2022,
  publisher    = {Zenodo},
  version      = {v0.1},
  doi          = {10.5281/zenodo.6415138},
  url          = {https://doi.org/10.5281/zenodo.6415138}
}
```

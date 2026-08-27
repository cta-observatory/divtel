# divtel: Divergent pointing mode for Imaging Atmospheric Cherenkov Telescopes arrays


<div align="left">

[![Build status](https://github.com/cta-observatory/divtel/actions/workflows/build.yml/badge.svg)](https://github.com/cta-observatory/divtel/actions/workflows/build.yml)
[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://pypi.org/project/divtel/)
[![Semantic Versions](https://img.shields.io/badge/%20%20%F0%9F%93%A6%F0%9F%9A%80-semantic--versions-e10079.svg)](https://github.com/cta-observatory/divtel/releases)
[![License](https://img.shields.io/github/license/cta-observatory/divtel?style=flat)](https://github.com/cta-observatory/divtel/blob/master/LICENSE)
[![Documentation](https://img.shields.io/github/actions/workflow/status/cta-observatory/divtel/doc.yml?branch=master&label=Documentation)](https://cta-observatory.github.io/divtel/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6415137.svg)](https://doi.org/10.5281/zenodo.6415137)
</div>

divtel makes toy simulations for the divergent pointing mode of Imaging Atmospheric Cherenkov Telescopes arrays.

Point an array's telescopes slightly away from one another and it sees a wider
patch of sky, but fewer telescopes see any given part of it, and a shower
needs at least two of them to be reconstructed stereoscopically. divtel lets
you set up that trade-off and measure both sides of it.

**[Try it in your browser](https://cta-observatory.github.io/divtel/examples.html)**,
sliders, no install.

## 👨‍💻 Install

```
pip install divtel
```

## 🚀 Quickstart

```python
import astropy.units as u
import matplotlib.pyplot as plt
from divtel.telescope import Telescope, Array
from divtel.visualization import display_hyper_fov

# Four telescopes on a 100 m square, each with a ~5.7 degree camera.
array = Array([
    Telescope(x * u.m, y * u.m, 0 * u.m, focal=20 * u.m, camera_radius=1 * u.m)
    for x, y in [(100, 0), (0, 100), (-100, 0), (0, -100)]
])

# Point them divergently around a mean direction of alt=70, az=180.
array.divergent_pointing(0.02, 70 * u.deg, 180 * u.deg)

# How much sky does the array see, and how much of it stereoscopically?
covered, patches = array.hyper_fov()        # 45.96 deg2
stereo, _ = array.hyper_fov(m_cut=2)        # 30.27 deg2, seen by 2+ telescopes

fig, (ground, sky) = plt.subplots(1, 2, figsize=(11, 5))
array.display_2d(projection="xy", ax=ground)
display_hyper_fov(array, ax=sky)
plt.show()
```

Pointed in parallel (`div=0`) the same array sees 25.7 deg2, all of it at
multiplicity 4. A `div` of 0.02 buys 79% more sky, of which 18% more is still
stereoscopic. But a shower now lands on two telescopes where it used to land
on four. That is the whole trade-off, and the
[user guide](https://cta-observatory.github.io/divtel/guide.html) walks through
it properly: the coordinate frame, what `div` really means, and how the hyper
field of view is computed.

## 🛠 Development

### Install from source

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

### Building the documentation

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

### Working on the interactive demo

`examples/marimo/interactive_display.py` is a [marimo](https://marimo.io)
notebook. Sphinx exports it to WebAssembly during the build, so the published
page ships its own Python interpreter and runs entirely in the reader's
browser — sliders included, with no server and nothing to install.

That export shells out to [uv](https://docs.astral.sh/uv/), which is why `uv`
is itself a documentation dependency.

The Jupyter version in `examples/notebooks/interactive_display.ipynb` is kept
for running locally, but is deliberately not built into the site: its
`ipywidgets` sliders need a Python kernel, so on a static page they would
render as controls that cannot move.


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

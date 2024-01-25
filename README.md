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

``` 
git clone -b irene_tests https://github.com/cta-observatory/divtel.git
cd divtel
conda env create  -f environment.yml
conda activate divtel
pip install -r requirements.txt
pip install -e ".[tests]"
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

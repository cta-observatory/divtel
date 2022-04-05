import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u

from .const import *

from .cta import CTA_Info

from .ArrayConfig import LoadConfig

from astropy.visualization import astropy_mpl_style, quantity_support
plt.style.use(astropy_mpl_style)

quantity_support()

__version__ = 0.1
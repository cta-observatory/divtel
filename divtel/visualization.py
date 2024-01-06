import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import widgets
from IPython.display import display

import copy
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroplan.plots import plot_sky
from astroplan import FixedTarget
from . import utils
from .const import COLORS
#from . import pointing

from matplotlib.transforms import Affine2D
from astropy.visualization.wcsaxes import SphericalCircle

import healpy as hp
import tqdm


def display_1d(table, proj, ax=None, labels=None, **kwargs):
    xb = utils.calc_mean(table, proj[0])
    
    ax = plt.figure().add_subplot(111) if ax is None else ax

    for i, [tels, label] in enumerate(zip(table.groups, labels)):
        c = COLORS(i)
        for val in tels[proj[0]]:
            ax.axvline(val, label=label, color=c, **kwargs)
            label='_nolegend_'

    ax.axvline(xb, color="r", label='barycenter', **kwargs)
    ax.set_xlabel(f"{proj[0]} [m]")
    ax.set_yticks([0, 1])
    ax.legend(frameon=False)

    return ax

def display_2d(table, proj, ax=None, labels=None, **kwargs):
    if ax is None:
        ax = plt.figure().add_subplot(111)
            
    scale = 1
    
    b_output = utils.calc_mean(table, [proj[0], proj[1], f"p_{proj[0]}", f"p_{proj[1]}"])
    
    for i, [tels, label] in enumerate(zip(table.groups, labels)):
        xx = tels[proj[0]]
        yy = tels[proj[1]]
        xv = tels[f"p_{proj[0]}"]
        yv = tels[f"p_{proj[1]}"]
        ids = tels["id"]

        s = ax.scatter(xx, yy, label=label, **kwargs)
        ax.quiver(xx, yy, xv, yv, color=s.get_facecolor())

        for i, x, y in zip(ids, xx, yy):
            ax.annotate(i, (x, y))
    
    ax.scatter(b_output[0], b_output[1], marker='+', label='barycenter', color="r")
    ax.quiver(*b_output, color="r")
    ax.set_xlabel(f"{proj[0]} [m]")
    ax.set_ylabel(f"{proj[1]} [m]")

    ax.grid('on')
    ax.axis('equal')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim(xlim[0] - 0.25 * np.abs(xlim[0]), xlim[1] + 0.25 * np.abs(xlim[1]))
    ax.set_ylim(ylim[0] - 0.25 * np.abs(ylim[0]), ylim[1] + 0.25 * np.abs(ylim[1]))
    ax.legend(frameon=False)

    return ax

def display_3d(table, proj, ax=None, labels=None, **kwargs):
    ax = plt.figure().add_subplot(111, projection='3d')

    scale = 1

    max_range = []
    for axis in ["x", "y", "z"]:
        max_range.append(table[axis].max() - table[axis].min())
    max_range = max(max_range)

    for i, [tels, label] in enumerate(zip(table.groups, labels)):
        xx = tels["x"]
        yy = tels["y"]
        zz = tels["z"]
        c = COLORS(i)
        ax.quiver(xx, yy, zz, 
                tels["p_x"], tels["p_y"], tels["p_z"],
                length=max_range,
                label=label,
                color=c,
                )

        Xb = scale * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + scale * (xx.max() + xx.min())
        Yb = scale * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + scale * (yy.max() + yy.min())
        Zb = scale * max_range * np.mgrid[-0.01:2:2, -0.01:2:2, -0.01:2:2][2].flatten() + scale * (zz.max() + zz.min())
        
        for xb, yb, zb in zip(Xb, Yb, Zb):
            ax.plot([xb], [yb], [zb], 'w')
        
    xx = utils.calc_mean(table, proj[0])
    yy = utils.calc_mean(table, proj[1])
    zz = utils.calc_mean(table, proj[2])
    xbv = utils.calc_mean(table, f"p_{proj[0]}")
    ybv = utils.calc_mean(table, f"p_{proj[1]}")
    zbv = utils.calc_mean(table, f"p_{proj[2]}")

    ax.quiver(xx, yy, zz, 
            xbv, ybv, zbv,
            color="r",
            length=max_range,
            label='barycenter',
            )
     
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    ax.legend(frameon=False)
    
    return ax

def display_barycenter(table, proj, ax=None, labels=None, fig=None, **kwargs):
    if fig is None:
        fig = plt.figure() 

    if ax is None:
        ax = fig.add_subplot(111) 

    scale = 1
    
    for i, (tab, label) in enumerate(zip(table.groups, labels)):
        output = utils.calc_mean(tab, [proj[0], proj[1], f"p_{proj[0]}", f"p_{proj[1]}"])
        s = ax.scatter(output[0], output[1], color=COLORS(i),)
        ax.quiver(*output, color=s.get_facecolor())

        ax.annotate(label, (output[0], output[1]))

    ax.set_xlabel(f"{proj[0]} [m]")
    ax.set_ylabel(f"{proj[1]} [m]")

    ax.grid('on')
    ax.axis('equal')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim(xlim[0] - 0.25 * np.abs(xlim[0]), xlim[1] + 0.25 * np.abs(xlim[1]))
    ax.set_ylim(ylim[0] - 0.25 * np.abs(ylim[0]), ylim[1] + 0.25 * np.abs(ylim[1]))
    ax.legend(frameon=False)

    return ax

def interactive_barycenter(array, proj="xy", overwrite=True, group=False):
    if overwrite:
        new_array = array
    else:
        new_array = copy.deepcopy(array)
    
    fig = plt.figure()

    def update(div=0, az=0, alt=0):
        new_array.divergent_pointing(div, az=az, alt=alt, units='deg')
        new_array.__convert_units__(toDeg=True)
        plt.cla()
        grouped_table, labels = new_array.group_by(group)
        display_barycenter(grouped_table, proj, labels=labels, fig=fig)
        fig.canvas.draw_idle()

    div_s = widgets.FloatLogSlider(value=new_array.div, base=10, min=-4, max=0, step=0.2, description='Divergence')
    az_s = widgets.FloatSlider(value=new_array.pointing["az"].value, min=0, max=360, step=0.01, description='Azimuth [deg]')
    alt_s = widgets.FloatSlider(value=new_array.pointing["alt"].value, min=0, max=90, step=0.01, description='Altitude [deg]')
    
    ui = widgets.HBox([div_s, alt_s, az_s])
    out = widgets.interactive_output(update, {'div': div_s, 'az': az_s, 'alt': alt_s})
    display(ui, out)

    return new_array


def multiplicity_plot(array, fig=None):
   
    nside = 512
    map_multiplicity = np.zeros(hp.nside2npix(nside), dtype=np.int8)
    counter = np.arange(0, hp.nside2npix(nside))
    ra, dec = hp.pix2ang(nside, counter, True, lonlat=True)
    coordinate = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    coord = array.get_pointing_coord(icrs=False)
    for i in tqdm.tqdm(range(len(array.telescopes))):
        pointing = SkyCoord(ra=coord.az[i].degree, dec=coord.alt[i].degree, unit='deg')
        r_fov = np.arctan((array.telescopes[i].camera_radius/array.telescopes[i].focal).to(u.dimensionless_unscaled)).to(u.deg)
        mask = coordinate.separation(pointing) < r_fov
        map_multiplicity[mask] += 1
    
    
    R=np.sqrt(array.hFoV()[0]/np.pi) + 5
    hp.cartview(map_multiplicity, rot=[array.pointing["az"].value, array.pointing["alt"].value],
                lonra=[-R,R], latra=[-R,R], nest=True, cmap='viridis', title=f"{array.frame.site} div={array.div}")
    # Annotate with axis labels:
    plt.annotate('Right Ascension (degrees)', xy=(0.5, -0.05), xycoords='axes fraction', ha='center', va='center')
    plt.annotate('Declination (degrees)', 
                 xy=(-0.05, 0.5), xycoords='axes fraction', 
                 ha='center', va='center', rotation='vertical')
    hp.graticule(dpar=5, dmer=5, coord='G', color='gray', lw=0.5)
    
    plt.show()


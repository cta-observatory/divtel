"""
Functions to define telescopes pointings
We use the same reference frame as simtel_array:
X is pointing North
Y is pointing West
Z is pointing upward
Az is taken clock-wise from X (towards Y) and between -180 and 180 degrees
Alt is taken from ground (towards Z) and between -90 and 90 degrees
"""

import numpy as np
import astropy.units as u


def alt_az_to_vector(alt, az):
    """
    Compute a pointing vector from an alt,az pointing direction

    Parameters
    ----------
    alt: angle in rad
    az: angle in rad

    Returns
    -------
    np.array([x, y, z])
    """
    x = np.cos(alt.to(u.rad)) * np.cos(az.to(u.rad))
    y = -np.cos(alt.to(u.rad)) * np.sin(az.to(u.rad))
    z = np.sin(alt.to(u.rad))
    return np.array([x, y, z])

def _norm_div(div, scale=100):
    """
    Transformation function from div parameter to norm to compute the position of g_point

     Parameters
    ----------
    div: float
    scale: float
        telescope distance from barycenter at which div = divergence_angle/90deg

    Returns
    -------
    float
    """
    return scale/np.tan(np.arcsin(div))

def pointG_position(barycenter, div, alt_mean, az_mean):
    """
    Compute the position of g_point for the pointing

    Parameters
    ----------
    barycenter: np.array([x,y,z])
        position of the barycenter of the array
    div: float
    alt_mean: `astropy.Quantity`
        mean pointing altitude in radians from which to diverge
    az_mean: `astropy.Quantity`
        mean pointing azimuth in radians from which to diverge

    Returns
    -------
    Numpy array [Gx, Gy, Gz]
    """
    norm = _norm_div(div)
    g_x = barycenter[0] - norm * np.cos(alt_mean) * np.cos(az_mean)
    g_y = barycenter[1] + norm * np.cos(alt_mean) * np.sin(az_mean)
    g_z = barycenter[2] - norm * np.sin(alt_mean)
    return np.array([g_x, g_y, g_z])


def tel_div_pointing(tel_position, g_point):
    """
    Divergent pointing to a point G.
    Update telescope pointing

    Parameters
    ----------
    tel_position: np.array([x, y, z])
        telescope coordinates
    g_point: numpy.array([Gx, Gy, Gz])
    """
    GT = np.sqrt(((tel_position - g_point) ** 2).sum())
    alt_tel = np.arcsin((tel_position[2] - g_point[2]) / GT)
    az_tel = np.arctan2(-(tel_position[1] - g_point[1]), (tel_position[0] - g_point[0]))
    return alt_tel, az_tel

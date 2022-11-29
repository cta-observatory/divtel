import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u

def calc_mean(table, columns):
    """
    Calculate a mean value of columns

    Parameters
    ----------
    table: astropy.table
    columns: str or array
        names of columns

    Returns
    -------
    float
    """

    if np.size(columns) == 1:
        columns = [columns]
    
    mean = []
    for column in columns:
        if column in ["p_x", "p_y", "p_z"]:
            mean.append(np.average(table[column], weights=table["d_tel"]))
        else:
            mean.append(np.mean(table[column]))
    return tuple(mean)
    
def deg2rad(table, toDeg=False):
    """
    Convert units in a table from rad to deg and vice versa.
    `az`, `alt`, `zn`, `radius`, `fov` units are changed.

    Parameters
    ----------
    table: astropy.table
    toDeg: bool, optional
        If True, the output table is in `deg`. Otherwise, in `rad`.

    Returns
    -------
    astropy.table
    """
    if toDeg:
        for par in ["az", "alt", "zn"]:
            table[par] = table[par].to(u.deg)
            table[par].info.format = '{:.3f}'

        table["radius"] = convert_radius(table["radius"], table["focal"], toDeg=toDeg)
        table["radius"].info.format = '{:.3f}'

        table["fov"]    = table["fov"].to(u.deg**2)
        table["fov"].info.format = '{:.3f}'
    else:
        for par in ["az", "alt", "zn"]:
            table[par] = table[par].to(u.rad)
            table[par].info.format = '{:.3f}'
        
        table["radius"] = convert_radius(table["radius"], table["focal"], toDeg=toDeg)
        table["radius"].info.format = '{:.3f}'
        
        table["fov"]    = table["fov"].to(u.rad**2)
        table["fov"].info.format = '{:.3f}'
    return table

def convert_radius(radius, focal, toDeg=False):
    """
    Convert the unit of camera radius from degree to meter, and vice versa.

    Parameters
    ----------
    radius: float, astropy.Quantity
        camera radius of a telecope
    focal: float, astropy.Quantity 
        focal length of a telescope
    toDeg: bool, optional
        If True, the output is in degree.
        If False, the output is in meter.

    Returns
    -------
    astropy.Quantity
    """
    if toDeg and radius.unit == u.m:
        temp = np.arctan(np.asarray(radius/focal))
        temp = temp*u.rad
        radius = temp.to(u.deg)
    elif not(toDeg) and radius.unit == u.deg:        
        temp = radius.to(u.rad)
        radius= np.tan(temp.value)*focal
    return radius


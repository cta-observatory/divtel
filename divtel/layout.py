"""
Read an array layout from file.

Layouts are ECSV tables -- plain text with a small YAML header naming each
column and its unit. Carrying the units in the file removes the ambiguity that
plain columns of numbers leave behind: a camera radius may be quoted as a length
or as an angle on the sky, and the two look identical once the unit is gone.
Here the unit decides, so there is no flag to get wrong.

A layout needs the columns ``x``, ``y``, ``z``, ``focal`` and ``camera_radius``,
and may carry an ``id`` column of telescope identifiers.
"""

import astropy.units as u
import numpy as np
from astropy.table import QTable

from .telescope import Array, Telescope

__all__ = ["load_array", "load_table", "camera_radius_in_metres"]

REQUIRED_COLUMNS = ("x", "y", "z", "focal", "camera_radius")


def camera_radius_in_metres(camera_radius, focal):
    """
    Camera radius as a length, whether it was given as one or as an angle.

    A camera of radius `r` at the focus of a telescope of focal length `f`
    subtends a half-angle `arctan(r/f)` on the sky, so an angular radius
    converts back with `r = f tan(angle)`. Angles are what instrument papers
    quote; metres are what the geometry needs.

    Parameters
    ----------
    camera_radius: `astropy.Quantity`
        radius of the camera, as a length or as an angle
    focal: `astropy.Quantity`
        focal length

    Returns
    -------
    `astropy.Quantity`
        radius as a length
    """
    physical_type = camera_radius.unit.physical_type

    if physical_type == "angle":
        return (focal * np.tan(camera_radius)).to(u.m)
    if physical_type == "length":
        return camera_radius.to(u.m)

    raise ValueError(
        f"camera_radius must be a length or an angle, got {camera_radius.unit} "
        f"({physical_type})"
    )


def load_table(path):
    """
    Read a layout file as a table, without building any telescope.

    Useful to look a layout over, or to edit it, before turning it into an
    `Array` -- `load_array` accepts the table this returns.

    Parameters
    ----------
    path: str or `pathlib.Path`

    Returns
    -------
    `astropy.table.QTable`
        with the columns of the file, as quantities

    Raises
    ------
    ValueError
        if a required column is missing or the table is empty
    """
    table = QTable.read(path, format="ascii.ecsv")

    missing = [name for name in REQUIRED_COLUMNS if name not in table.colnames]
    if missing:
        raise ValueError(
            f"{path} is missing the column(s) {', '.join(missing)}; a layout needs "
            f"{', '.join(REQUIRED_COLUMNS)}"
        )

    if len(table) == 0:
        raise ValueError(f"{path} holds no telescope")

    return table


def load_array(path):
    """
    Build an `Array` from a layout file.

    Parameters
    ----------
    path: str, `pathlib.Path` or `astropy.table.QTable`
        an ECSV layout file, or a table already read from one

    Returns
    -------
    `Array`

    Notes
    -----
    Telescope ids come from the file's ``id`` column when it has one, so that a
    layout read twice gives the same ids both times and those ids mean what the
    layout says they mean. Without that column telescopes are numbered 1 to N in
    file order.

    Examples
    --------
    >>> from importlib.resources import files
    >>> from divtel.layout import load_array
    >>> array = load_array(files("divtel") / "data" / "cta-north-lapalma-alpha-prod6.ecsv")
    >>> len(array.telescopes)
    13
    """
    table = path if isinstance(path, QTable) else load_table(path)

    if "id" in table.colnames:
        ids = [int(tel_id) for tel_id in table["id"]]
        duplicates = sorted({tel_id for tel_id in ids if ids.count(tel_id) > 1})
        if duplicates:
            raise ValueError(
                f"{path} reuses the telescope id(s) {', '.join(map(str, duplicates))}; "
                "ids must be unique"
            )
    else:
        ids = range(1, len(table) + 1)

    telescopes = [
        Telescope(
            row["x"],
            row["y"],
            row["z"],
            row["focal"],
            camera_radius_in_metres(row["camera_radius"], row["focal"]),
            tel_id=tel_id,
        )
        for row, tel_id in zip(table, ids, strict=True)
    ]

    return Array(telescopes)

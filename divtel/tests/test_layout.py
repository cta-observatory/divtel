"""Layouts are read from ECSV files that carry their own units."""

from importlib.resources import files

import astropy.units as u
import pytest

from divtel.layout import camera_radius_in_metres, load_array, load_table

DATA = files("divtel") / "data"

HEADER = """# %ECSV 1.0
# ---
# datatype:
# - {{name: x, datatype: float64, unit: m}}
# - {{name: y, datatype: float64, unit: m}}
# - {{name: z, datatype: float64, unit: m}}
# - {{name: focal, datatype: float64, unit: m}}
# - {{name: camera_radius, datatype: float64, unit: {radius_unit}}}
x y z focal camera_radius
"""


def write_layout(tmp_path, body, radius_unit="m", header=HEADER):
    path = tmp_path / "layout.ecsv"
    path.write_text(header.format(radius_unit=radius_unit) + body)
    return path


def test_the_packaged_dummy_array_loads():
    array = load_array(DATA / "dummy_array.ecsv")

    assert len(array.telescopes) == 5
    assert array.telescopes[0].focal == 28 * u.m


def test_the_packaged_la_palma_layout_loads():
    array = load_array(DATA / "la_palma_4LST_15MST.ecsv")

    assert len(array.telescopes) == 19
    # four LSTs of 28 m focal, then fifteen MSTs of 16 m
    focals = [tel.focal for tel in array.telescopes]
    assert focals[:4] == [28 * u.m] * 4
    assert focals[4:] == [16 * u.m] * 15


def test_positions_are_quantities():
    telescope = load_array(DATA / "la_palma_4LST_15MST.ecsv").telescopes[0]

    assert isinstance(telescope.position, u.Quantity)
    assert telescope.position.unit.is_equivalent(u.m)


def test_ids_come_from_the_file():
    """CTAO ids carry meaning -- LSTs 1-4, MSTs 5 on -- so the file sets them."""
    array = load_array(DATA / "la_palma_4LST_15MST.ecsv")

    assert [tel.id for tel in array.telescopes] == list(range(1, 20))


def test_ids_are_the_same_every_time_a_layout_is_loaded():
    """Regression test for issue #13: ids used to come from a global counter."""
    first = load_array(DATA / "la_palma_4LST_15MST.ecsv")
    second = load_array(DATA / "la_palma_4LST_15MST.ecsv")

    assert [tel.id for tel in first.telescopes] == [tel.id for tel in second.telescopes]


def test_a_layout_without_ids_is_numbered_in_file_order(tmp_path):
    path = write_layout(tmp_path, "0 0 0 28 1\n10 0 0 28 1\n20 0 0 28 1\n")

    assert [tel.id for tel in load_array(path).telescopes] == [1, 2, 3]


def test_duplicate_ids_are_rejected(tmp_path):
    header = HEADER.replace(
        "# datatype:", "# datatype:\n# - {{name: id, datatype: int64}}"
    ).replace("x y z focal", "id x y z focal")
    path = write_layout(tmp_path, "1 0 0 0 28 1\n1 10 0 0 28 1\n", header=header)

    with pytest.raises(ValueError, match="reuses the telescope id"):
        load_array(path)


def test_an_angular_camera_radius_is_converted_to_a_length():
    """A camera radius may be quoted as a length or as an angle; the unit decides."""
    in_metres = load_array(DATA / "la_palma_4LST_15MST.ecsv")
    in_degrees = load_array(DATA / "la_palma_4LST_15MST_deg.ecsv")

    for metres, degrees in zip(in_metres.telescopes, in_degrees.telescopes, strict=True):
        assert degrees.camera_radius.unit.is_equivalent(u.m)
        assert abs(metres.camera_radius - degrees.camera_radius) < 0.1 * u.mm


def test_camera_radius_conversion_round_trips():
    radius = camera_radius_in_metres(2.1476 * u.deg, 28 * u.m)

    assert radius.unit.is_equivalent(u.m)
    assert radius.to_value(u.m) == pytest.approx(1.05, abs=1e-4)


def test_a_camera_radius_that_is_neither_a_length_nor_an_angle_is_rejected():
    with pytest.raises(ValueError, match="must be a length or an angle"):
        camera_radius_in_metres(1.0 * u.s, 28 * u.m)


def test_a_missing_column_is_named(tmp_path):
    header = HEADER.replace(
        "# - {{name: z, datatype: float64, unit: m}}\n", ""
    ).replace("x y z focal", "x y focal")
    path = write_layout(tmp_path, "0 0 28 1\n", header=header)

    with pytest.raises(ValueError, match="missing the column"):
        load_table(path)


def test_an_empty_layout_is_rejected(tmp_path):
    path = write_layout(tmp_path, "")

    with pytest.raises(ValueError, match="holds no telescope"):
        load_table(path)


def test_a_table_can_be_read_edited_and_then_built():
    table = load_table(DATA / "la_palma_4LST_15MST.ecsv")

    assert len(load_array(table[:4]).telescopes) == 4

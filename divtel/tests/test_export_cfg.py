import astropy.units as u
import pytest

from divtel.telescope import Array, Telescope


@pytest.fixture
def square_array():
    """Four telescopes on a 100 m square."""
    return Array(
        [
            Telescope(100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            Telescope(0 * u.m, 100 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            Telescope(-100 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
            Telescope(0 * u.m, -100 * u.m, 0 * u.m, 20 * u.m, 1 * u.m),
        ]
    )


def test_export_cfg_needs_a_pointing(square_array, tmp_path):
    """Exporting an array that was never pointed is a mistake, not a 0/0 file."""
    with pytest.raises(ValueError, match="no pointing yet"):
        square_array.export_cfg(outdir=tmp_path)


def test_export_cfg_writes_one_block_per_telescope(square_array, tmp_path):
    square_array.divergent_pointing(0.05, 70 * u.deg, 180 * u.deg)
    path = square_array.export_cfg(outdir=tmp_path)

    text = path.read_text()
    assert path.exists()
    for n in range(1, 5):
        assert f"#elif TELESCOPE == {n}\n" in text
    assert "#elif TELESCOPE == 5" not in text
    assert text.count("TELESCOPE_THETA") == 5  # block 0 plus one per telescope
    assert text.rstrip().endswith("array_trigger = array_trigger_ultra6_diver-test.dat")


def test_export_cfg_default_filename_records_the_pointing(square_array, tmp_path):
    square_array.divergent_pointing(0.05, 70 * u.deg, 180 * u.deg)
    path = square_array.export_cfg(outdir=tmp_path)
    assert path.name == "CTA-ULTRA6-LaPalma-div0_05-az180_0-alt70_0.cfg"


def test_export_cfg_converts_to_simtel_angle_convention(square_array, tmp_path):
    """theta = 90 - alt, phi = (360 - az) mod 360, with phi never negative."""
    square_array.divergent_pointing(0, 70 * u.deg, 180 * u.deg)
    text = square_array.export_cfg(outdir=tmp_path).read_text()

    # Parallel pointing: every telescope sits at the mean pointing.
    assert text.count("TELESCOPE_THETA=20.00") == 5
    assert text.count("TELESCOPE_PHI=180.00") == 5


def test_export_cfg_keeps_phi_in_zero_to_360(square_array, tmp_path):
    """divtel azimuths run -180..180; the flipped angle must be wrapped."""
    square_array.divergent_pointing(0.2, 70 * u.deg, 10 * u.deg)
    text = square_array.export_cfg(outdir=tmp_path).read_text()

    phis = [float(line.split("=")[1])
            for line in text.splitlines() if "TELESCOPE_PHI" in line]
    assert all(0 <= phi < 360 for phi in phis), phis


def test_export_cfg_defaults_to_four_lsts_then_msts(square_array, tmp_path):
    square_array.divergent_pointing(0, 70 * u.deg, 180 * u.deg)
    text = square_array.export_cfg(outdir=tmp_path).read_text()
    assert text.count("#  include <CTA-ULTRA6-LST.cfg>") == 5  # block 0 + 4 LSTs
    assert "NectarCam" not in text


def test_export_cfg_accepts_explicit_tel_configs(square_array, tmp_path):
    square_array.divergent_pointing(0, 70 * u.deg, 180 * u.deg)
    configs = ["A.cfg", "B.cfg", "C.cfg", "D.cfg"]
    text = square_array.export_cfg(outdir=tmp_path, tel_configs=configs).read_text()
    for cfg in configs:
        assert f"#  include <{cfg}>" in text


def test_export_cfg_rejects_mismatched_tel_configs(square_array, tmp_path):
    square_array.divergent_pointing(0, 70 * u.deg, 180 * u.deg)
    with pytest.raises(ValueError, match="2 entries but the array has 4"):
        square_array.export_cfg(outdir=tmp_path, tel_configs=["A.cfg", "B.cfg"])

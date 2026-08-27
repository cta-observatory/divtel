import divtel

def test_version():
    divtel.__version__


def test_fov_convertible_to_steradian():
    """`fov` is a solid angle, so it must convert to `u.sr`."""
    import astropy.units as u
    from divtel.telescope import Telescope

    telescope = Telescope(0 * u.m, 0 * u.m, 0 * u.m, 28 * u.m, 1 * u.m)
    telescope.fov.to(u.sr)


def test_dummy_array_data_loads():
    """The packaged `data/dummy_array.txt` loads into a working Array."""
    import numpy as np
    import astropy.units as u
    from importlib.resources import files
    from divtel.telescope import Array, Telescope

    path = files("divtel") / "data" / "dummy_array.txt"
    rows = np.genfromtxt(path, delimiter=",")

    array = Array([
        Telescope(x * u.m, y * u.m, z * u.m, focal * u.m, camera_radius * u.m)
        for x, y, z, focal, camera_radius in rows
    ])

    assert len(array.telescopes) == 5


def test_point_to_object_round_trip():
    """`pointing_vector` must point at the object it was aimed at."""
    import numpy as np
    import astropy.units as u
    from divtel.telescope import Telescope

    telescope = Telescope(0 * u.m, 0 * u.m, 0 * u.m, 28 * u.m, 1 * u.m)
    obj = np.array([1000., 1000., 1000.])
    telescope.point_to_object(obj)

    np.testing.assert_allclose(telescope.pointing_vector,
                               obj / np.linalg.norm(obj),
                               atol=1e-12)


def test_point_to_object_round_trip_off_origin():
    """The round trip must hold for a telescope away from the origin too."""
    import numpy as np
    import astropy.units as u
    from divtel.telescope import Telescope

    position = np.array([-40., 75., 12.])
    telescope = Telescope(*(position * u.m), 28 * u.m, 1 * u.m)

    for obj in (np.array([300., -200., 900.]),
                np.array([-500., 640., 1200.]),
                np.array([120., 0., 400.])):
        telescope.point_to_object(obj)
        direction = obj - position
        np.testing.assert_allclose(telescope.pointing_vector,
                                   direction / np.linalg.norm(direction),
                                   atol=1e-12)

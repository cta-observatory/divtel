==========
User guide
==========

divtel is a small toy model of an Imaging Atmospheric Cherenkov Telescope
(IACT) array. You place telescopes on the ground, point them -- in parallel or
divergently -- and measure the sky they collectively see.

This guide walks through that from the start. For the signatures and every
argument, see the :doc:`API documentation <docstring>`; for a version you can
play with in the browser, see :doc:`examples`.

Installation
============

.. code-block:: console

    pip install divtel

The examples below also use ``matplotlib``, which comes with divtel.


The coordinate frame
====================

divtel uses the same ground frame as ``simtel_array``:

* **x** points North, **y** points West, **z** points up. Positions are
  lengths, normally metres.
* **azimuth** is measured clockwise from x towards y, in ``[-180, 180]``
  degrees.
* **altitude** is measured from the ground towards z, in ``[-90, 90]``
  degrees. Altitude 90 is the zenith.

Every quantity that carries a physical dimension is an
`astropy.units.Quantity <https://docs.astropy.org/en/stable/units/>`_. This is
not decoration: divtel mixes metres and angles freely, and the units are what
keeps that honest. Pass ``100 * u.m``, not ``100``.


Building an array
=================

A :class:`~divtel.telescope.Telescope` is a ground position plus the two
numbers that set how much sky it sees: the focal length and the camera radius.
An :class:`~divtel.telescope.Array` is a list of them.

.. code-block:: python

    import astropy.units as u
    from divtel.telescope import Telescope, Array

    # Four telescopes on a 100 m square -- roughly the H.E.S.S.-I layout.
    array = Array([
        Telescope(x * u.m, y * u.m, 0 * u.m,
                  focal=20 * u.m, camera_radius=1 * u.m)
        for x, y in [(100, 0), (0, 100), (-100, 0), (0, -100)]
    ])

The two camera numbers only ever enter as their ratio, which fixes the angular
radius of the disc a telescope sees on the sky:

.. code-block:: python

    >>> array.telescopes[0].fov_radius.to(u.deg)
    <Quantity 2.86240523 deg>
    >>> array.telescopes[0].fov.to(u.deg**2)
    <Quantity 25.78310078 deg2>

``fov_radius`` is the half-angle subtended by the camera; ``fov`` is the solid
angle of the whole disc. A 20 m focal length with a 1 m camera radius gives a
camera a little under 6 degrees across.

An array knows a few things about itself as a whole:

.. code-block:: python

    >>> array.barycenter                 # mean telescope position
    <Quantity [0., 0., 0.] m>
    >>> array.positions_array.shape      # one [x, y, z] row per telescope
    (4, 3)


Loading an array from file
==========================

Typing telescopes out gets old quickly. :func:`~divtel.layout.load_array` reads
a layout from an ECSV file -- plain text with a small header naming each column
and its unit:

.. code-block:: python

    from importlib.resources import files
    from divtel.layout import load_array

    array = load_array(files("divtel") / "data" / "la_palma_4LST_15MST.ecsv")

Two layouts ship with the package: ``dummy_array.ecsv``, five telescopes to try
things on, and ``la_palma_4LST_15MST.ecsv``, a CTA La Palma layout of 4 LSTs and
15 MSTs. A layout file needs the columns ``x``, ``y``, ``z``, ``focal`` and
``camera_radius``, and may carry an ``id`` column:

.. code-block:: text

    # %ECSV 1.0
    # ---
    # datatype:
    # - {name: id, datatype: int64}
    # - {name: x, datatype: float64, unit: m}
    # - {name: y, datatype: float64, unit: m}
    # - {name: z, datatype: float64, unit: m}
    # - {name: focal, datatype: float64, unit: m}
    # - {name: camera_radius, datatype: float64, unit: m}
    id x y z focal camera_radius
    1 -70.04 -7.23 54 28 1.05

Because the units live in the file, ``camera_radius`` may be given either as a
length or as the half-angle the camera subtends on the sky -- whichever the
instrument is quoted in. An angle is converted with ``radius = focal * tan(angle)``
on the way in, so the two spellings describe the same telescope:

.. code-block:: python

    >>> from divtel.layout import camera_radius_in_metres
    >>> camera_radius_in_metres(2.1476 * u.deg, 28 * u.m)
    <Quantity 1.05000713 m>

The ``id`` column sets :attr:`~divtel.telescope.Telescope.id`. CTAO numbers
telescopes by type -- LSTs 1 to 4, MSTs 5 to 14, telescopes added to try out a
configuration from 15 on -- so ids carry meaning, and
:meth:`~divtel.telescope.Array.group_by` selects on them. Without the column,
telescopes are numbered 1 to N in file order. Telescopes built by hand instead
take ids from a counter shared across the session, which makes them unique but
not meaningful; pass ``tel_id`` to :class:`~divtel.telescope.Telescope` to set
them yourself.

To look a layout over, or edit it, before building anything, read it as a table
first -- :func:`~divtel.layout.load_array` takes one:

.. code-block:: python

    from divtel.layout import load_table

    table = load_table(files("divtel") / "data" / "la_palma_4LST_15MST.ecsv")
    lsts = load_array(table[:4])


Pointing the array
==================

Parallel pointing
-----------------

Every telescope aimed the same way. In divtel this is just divergent pointing
with the divergence turned off:

.. code-block:: python

    array.divergent_pointing(0, 70 * u.deg, 180 * u.deg)

All four cameras now cover the same disc of sky: the array sees 25.7 deg²,
and every part of it is seen by all four telescopes.

Divergent pointing
------------------

Turn the knob up and the telescopes fan out:

.. code-block:: python

    array.divergent_pointing(0.02, 70 * u.deg, 180 * u.deg)

.. code-block:: python

    >>> array.pointing_altaz.to(u.deg)
    <Quantity [[  71.08431390,  180.        ],
               [  69.96853681, -176.65271002],
               [  68.93041909,  180.        ],
               [  69.96853681,  176.65271002]] deg>

The mean pointing is still (70°, 180°), but the four telescopes are now spread
about two degrees apart on the sky.

.. _the-div-parameter:

What ``div`` actually means
---------------------------

``div`` runs from 0 (parallel) to 1, and it works through a point **G** placed
*behind* the array. Every telescope is aimed directly away from G, so the
further away G sits, the more nearly parallel the array is:

.. code-block:: python

    from divtel import pointing

    >>> pointing.pointG_position(array.barycenter, 0.02,
    ...                          70 * u.deg, 180 * u.deg).to(u.km)
    <Quantity [ 1.70975866e+00,  2.09385047e-16, -4.69752332e+00] km>

At ``div = 0`` G recedes to infinity and the telescopes point in parallel; at
``div = 1`` it sits in the middle of the array and they point straight outward.

The one thing to keep in mind is that ``div`` is an engineering knob, not an
angle. It is defined through a fixed 100 m reference distance: ``div`` is the
sine of the divergence angle picked up by a telescope 100 m from the
barycenter, *measured perpendicular to the pointing direction*. So the same
``div`` produces a different real divergence on an array of a different size or
layout -- and even on the same array pointed a different way. If you need the
angle the array actually achieved, read it off ``pointing_altaz`` rather than
inferring it from ``div``.

Values outside ``[0, 1]`` raise a ``ValueError``; the geometry stops meaning
anything there.

The array remembers what it was last asked for, which saves carrying the
numbers around yourself:

.. code-block:: python

    >>> array.div
    0.02
    >>> array.mean_pointing
    (<Quantity 70. deg>, <Quantity 180. deg>)

``mean_pointing`` raises a ``ValueError`` on an array that has never been
pointed; ``div`` reports 0, which is what a fresh array's telescopes are
in fact doing -- they all sit at alt = az = 0 until you point them.

Pointing at a target
--------------------

A single telescope can also be aimed at a position in the ground frame, which
is useful for pointing at something at a finite distance:

.. code-block:: python

    array.telescopes[0].point_to_object([0, 0, 10000] * u.m)

Or directly in alt/az:

.. code-block:: python

    array.telescopes[0].point_to_altaz(70 * u.deg, 180 * u.deg)


The hyper field of view
=======================

This is the point of the whole exercise. Each telescope sees a disc on the
sky. Pointing divergently spreads those discs out, which buys **width** at the
cost of **depth**: a shower needs to be seen by at least two telescopes to be
reconstructed stereoscopically, and divergence trades away exactly that
overlap.

:meth:`~divtel.telescope.Array.hyper_fov` measures both sides of the bargain.
It cuts the union of the discs into patches and labels each with its
**multiplicity** -- how many telescopes see it:

.. code-block:: python

    >>> array.divergent_pointing(0.02, 70 * u.deg, 180 * u.deg)
    >>> total, patches = array.hyper_fov()          # everything covered
    >>> stereo, _ = array.hyper_fov(m_cut=2)        # seen by 2 or more
    >>> total, stereo
    (<Quantity 45.95557740 deg2>, <Quantity 30.27381956 deg2>)

Compare that with the parallel case, where all four discs coincide: 25.7 deg²
covered, all of it at multiplicity 4. A very modest ``div`` of 0.02 has bought
79% more sky, and 18% more of it is still stereoscopic -- but a shower now
lands on two telescopes where it used to land on four.

The second return value is the list of ``(polygon, multiplicity)`` patches, if
you want to do your own accounting:

.. code-block:: python

    >>> for polygon, multiplicity in patches:
    ...     print(multiplicity, round(polygon.area, 2))

Coordinates are degrees of offset from the array's mean pointing, in a Lambert
azimuthal equal-area projection centred on it. That projection is chosen so
patch areas are true solid angles whatever the array is pointing at, and so
that nothing breaks near the zenith -- where telescopes a fraction of a degree
apart on the sky are hundreds of degrees apart in azimuth.

Because the cameras here are about 5.7 degrees across while divergence swings
the pointings by tens of degrees, the discs come apart within the first few
hundredths of the slider. Push ``div`` past roughly 0.08 on this array and it
covers four separate patches of sky with no stereoscopic overlap at all.


Plotting
========

Two views, which together tell the story:

.. code-block:: python

    import matplotlib.pyplot as plt
    from divtel.visualization import display_hyper_fov

    fig, (ground, sky) = plt.subplots(1, 2, figsize=(11, 5))
    array.display_2d(projection="xy", ax=ground)
    display_hyper_fov(array, ax=sky)
    plt.show()

:meth:`~divtel.telescope.Array.display_2d` draws the telescopes on the ground
with an arrow each for where it points; pass ``projection='xz'`` or ``'yz'``
for a side view. :func:`~divtel.visualization.display_hyper_fov` draws the sky,
shaded by multiplicity, with the covered area in the title. Pass ``m_cut=2``
to have that title report the stereoscopic area instead -- the patches below
the cut are still drawn, faded, so the shape of the whole coverage stays
visible.

Both take an existing ``ax`` and return it, so they compose with whatever else
you are plotting.


Exporting to sim_telarray
=========================

Once a pointing looks right, :meth:`~divtel.telescope.Array.export_cfg` writes
it out as a ``sim_telarray`` configuration file, so the layout can be handed to
a real shower simulation:

.. code-block:: python

    >>> array.export_cfg(outdir="configs/", tel_configs=["lst.cfg"] * 4)
    PosixPath('configs/CTA-ULTRA6-LaPalma-div0_02-az180_0-alt70_0.cfg')

The default filename records the divergence and mean pointing, so a scan over
``div`` lands in distinct files without naming them. ``tel_configs`` gives the
per-telescope config each block ``#include``\ s, in array order; left out, it
defaults to the La Palma layout the writer was built for -- four LSTs followed
by MSTs with NectarCam -- which is only right if your array is that one.

Note that the angle conventions differ from divtel's: ``sim_telarray`` takes a
zenith angle rather than an altitude, and its azimuth runs the other way round.
The conversion is done for you on the way out.


Where to go next
================

* :doc:`examples` -- the same model with sliders, running in your browser.
* :doc:`docstring` -- the full API.
* `Source and issue tracker <https://github.com/cta-observatory/divtel>`_.

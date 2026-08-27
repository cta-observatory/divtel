==========
User guide
==========

divtel is a toy model of an Imaging Atmospheric Cherenkov Telescope (IACT)
array. Place telescopes on the ground, point them in parallel or
divergently, and measure the sky they see.

For full signatures, see the :doc:`API docs <docstring>`. For something to
play with in the browser, see :doc:`examples`.

Installation
============

divtel is on `PyPI <https://pypi.org/project/divtel/>`_:

.. code-block:: console

    pip install divtel

The examples below also use ``matplotlib``, which comes with divtel.


The coordinate frame
====================

divtel uses the same ground frame as ``simtel_array``:

* **x** points North, **y** points West, **z** points up, in metres.
* **azimuth** is measured clockwise from x towards y, ``[-180, 180]`` degrees.
* **altitude** is measured from the ground up, ``[-90, 90]`` degrees; 90 is
  the zenith.

Every physical quantity is an `astropy.units.Quantity
<https://docs.astropy.org/en/stable/units/>`_. divtel mixes metres and angles
freely, so pass ``100 * u.m``, not ``100``.


Building an array
=================

A :class:`~divtel.telescope.Telescope` is a ground position plus a focal
length and camera radius, which set how much sky it sees. An
:class:`~divtel.telescope.Array` is a list of them.

.. code-block:: python

    import astropy.units as u
    from divtel.telescope import Telescope, Array

    # Four telescopes on a 100 m square, roughly the H.E.S.S.-I layout.
    array = Array([
        Telescope(x * u.m, y * u.m, 0 * u.m,
                  focal=20 * u.m, camera_radius=1 * u.m)
        for x, y in [(100, 0), (0, 100), (-100, 0), (0, -100)]
    ])

The two camera numbers only ever enter as their ratio, the angular radius of
the disc a telescope sees on the sky:

.. code-block:: python

    >>> array.telescopes[0].fov_radius.to(u.deg)
    <Quantity 2.86240523 deg>
    >>> array.telescopes[0].fov.to(u.deg**2)
    <Quantity 25.78310078 deg2>

``fov_radius`` is the half-angle subtended by the camera; ``fov`` is the solid
angle of the whole disc. A 20 m focal length with a 1 m camera radius gives a
camera a little under 6 degrees across.

An array also knows its own geometry:

.. code-block:: python

    >>> array.barycenter                 # mean telescope position
    <Quantity [0., 0., 0.] m>
    >>> array.positions_array.shape      # one [x, y, z] row per telescope
    (4, 3)


Loading an array from file
==========================

Typing telescopes out by hand gets old fast.
:func:`~divtel.layout.load_array` reads a layout from an ECSV file instead,
plain text with a small header naming each column and its unit:

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

``camera_radius`` may be given as a length or as the half-angle the camera
subtends on the sky, since the units live in the file.
``radius = focal * tan(angle)`` converts the angle on the way in, so both
spellings describe the same telescope:

.. code-block:: python

    >>> from divtel.layout import camera_radius_in_metres
    >>> camera_radius_in_metres(2.1476 * u.deg, 28 * u.m)
    <Quantity 1.05000713 m>

The ``id`` column sets :attr:`~divtel.telescope.Telescope.id`. CTAO ids are
semantic: LSTs are 1-4, MSTs 5-14, and anything added to try out a
configuration is 15+. :meth:`~divtel.telescope.Array.group_by` selects on
these ids, so keep them meaningful. Without an ``id`` column, telescopes are
numbered 1 to N in file order; telescopes built by hand get ids from a
session-wide counter, unique but not meaningful, unless you pass ``tel_id``
yourself.

To inspect or edit a layout before building anything, read it as a table:
:func:`~divtel.layout.load_array` also accepts one:

.. code-block:: python

    from divtel.layout import load_table

    table = load_table(files("divtel") / "data" / "la_palma_4LST_15MST.ecsv")
    lsts = load_array(table[:4])


Pointing the array
==================

Parallel pointing
-----------------

Every telescope aimed the same way, divergence turned off:

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

``div`` runs from 0 (parallel) to 1. Internally it works through a point
**G** placed *behind* the array: every telescope aims directly away from G,
so the farther G sits, the more parallel the array is.

.. code-block:: python

    from divtel import pointing

    >>> pointing.pointG_position(array.barycenter, 0.02,
    ...                          70 * u.deg, 180 * u.deg).to(u.km)
    <Quantity [ 1.70975866e+00,  2.09385047e-16, -4.69752332e+00] km>

At ``div = 0`` G recedes to infinity and the telescopes point in parallel; at
``div = 1`` it sits in the middle of the array and they point straight outward.

``div`` is an engineering knob, not an angle: it's the sine of the
divergence picked up by a telescope 100 m from the barycenter, measured
perpendicular to the pointing direction. The same ``div`` value gives a
different real divergence on a different array, or even the same array
pointed elsewhere. Read ``pointing_altaz`` for the angle actually achieved.

Values outside ``[0, 1]`` raise a ``ValueError``; the geometry stops meaning
anything there.

Divergence spreads telescopes on both sides of the mean pointing, so pointing
near the horizon can push some of them underground. A flat array pointed at
the horizon won't show this, since there's no slope to tip telescopes past
it, but a sloped one will:

.. code-block:: python

    >>> hillside = Array([
    ...     Telescope(100 * u.m, 0 * u.m,  50 * u.m, 20 * u.m, 1 * u.m),
    ...     Telescope(-100 * u.m, 0 * u.m, -50 * u.m, 20 * u.m, 1 * u.m),
    ... ])
    >>> hillside.divergent_pointing(0.3, 5 * u.deg, 180 * u.deg)
    >>> np.round(hillside.pointing_altaz[:, 0].to(u.deg), 2)
    <Quantity [19.72,  0.  ] deg>

Telescope one, up the slope, ends up at 19.72°. Telescope two, down the
slope, would reach -3.06°, confirmed directly against
``pointing.tel_div_pointing``, which doesn't clamp. A real mount can't follow
the ground below the horizon, though, so ``divergent_pointing`` clamps it
there instead, at the azimuth it was already using. ``mean_pointing`` still
reports the direction you asked to diverge from, not what any one telescope
achieved, clamped or not.

The array remembers its last pointing:

.. code-block:: python

    >>> array.div
    0.02
    >>> array.mean_pointing
    (<Quantity 70. deg>, <Quantity 180. deg>)

``mean_pointing`` raises ``ValueError`` before the array has ever been
pointed. ``div`` reports 0 by default, since a fresh array's telescopes all
sit at alt = az = 0.

Pointing at a target
--------------------

A single telescope can also aim at a position in the ground frame, useful for
a source at finite distance:

.. code-block:: python

    array.telescopes[0].point_to_object([0, 0, 10000] * u.m)

Or directly in alt/az:

.. code-block:: python

    array.telescopes[0].point_to_altaz(70 * u.deg, 180 * u.deg)


The hyper field of view
=======================

Each telescope sees a disc on the sky. Divergent pointing spreads those discs
apart, trading **depth** for **width**: a shower needs at least two
telescopes to be reconstructed stereoscopically, and divergence trades away
that overlap.

:meth:`~divtel.telescope.Array.hyper_fov` measures both sides of the trade.
It cuts the union of the discs into patches and labels each with its
**multiplicity**, how many telescopes see it:

.. code-block:: python

    >>> array.divergent_pointing(0.02, 70 * u.deg, 180 * u.deg)
    >>> total, patches = array.hyper_fov()          # everything covered
    >>> stereo, _ = array.hyper_fov(m_cut=2)        # seen by 2 or more
    >>> total, stereo
    (<Quantity 45.95557740 deg2>, <Quantity 30.27381956 deg2>)

Compare that to parallel pointing, where all four discs coincide: 25.7 deg²
covered, all of it at multiplicity 4. A ``div`` of just 0.02 buys 79% more
sky, 18% of it still stereoscopic, but a shower now lands on two telescopes
where it used to land on four.

The second return value is the list of ``(polygon, multiplicity)`` patches, if
you want to do your own accounting:

.. code-block:: python

    >>> for polygon, multiplicity in patches:
    ...     print(multiplicity, round(polygon.area, 2))

Coordinates are degrees of offset from the array's mean pointing, in a
Lambert azimuthal equal-area projection: x is azimuth offset, y is altitude
offset. Patch areas come out as true solid angles at any pointing, including
near the zenith, where telescopes a fraction of a degree apart on the sky can
be hundreds of degrees apart in azimuth.

Because the cameras here are about 5.7 degrees across while divergence swings
the pointings by tens of degrees, the discs come apart within the first few
hundredths of the slider. Push ``div`` past roughly 0.08 on this array and it
covers four separate patches of sky with no stereoscopic overlap at all.


How well the sky is covered
===========================

The hyper FoV says how much sky is covered.
:meth:`~divtel.telescope.Array.multiplicity_profile` says how well, how much
of it is seen by exactly one telescope, by exactly two, and so on:

.. code-block:: python

    >>> array = load_array(files("divtel") / "data" / "la_palma_4LST_15MST.ecsv")
    >>> array.divergent_pointing(0.05, 70 * u.deg, 180 * u.deg)
    >>> multiplicity, area = array.multiplicity_profile()
    >>> multiplicity
    array([1, 2, 3, 4, 5])
    >>> np.round(area, 1)
    <Quantity [242. , 131.3,  70.8,   8. ,   0.6] deg2>

:meth:`~divtel.telescope.Array.multiplicity_moments` boils that down to two
numbers:

.. code-block:: python

    >>> mean, variance = array.multiplicity_moments()
    >>> f"{mean:.2f} +- {np.sqrt(variance):.2f}"
    '1.66 +- 0.81'


How the mean is defined
-----------------------

:meth:`~divtel.telescope.Array.hyper_fov` cuts the sky into disjoint patches,
each with a well-defined multiplicity :math:`m_i` and area :math:`A_i`. Mean
and variance are those values **weighted by solid angle**:

.. math::

   \langle m \rangle = \frac{\sum_i m_i A_i}{\sum_i A_i}
   \qquad
   \mathrm{Var}(m) = \frac{\sum_i \left(m_i - \langle m \rangle\right)^2 A_i}{\sum_i A_i}

Weighting by area matters: overlapping discs cut out slivers as often as
large regions, and an unweighted mean would count a sliver the same as the
whole field. ``hyper_fov`` works in an equal-area projection, so the areas
are true solid angles, and the patches partition the covered sky exactly, so
nothing is double-counted.

Working the example above by hand:

.. code-block:: text

    m :  [1       2       3      4     5   ]
    A :  [242.02  131.34  70.81  7.99  0.55]  deg²    sum = 452.70

    (1*242.02 + 2*131.34 + 3*70.81 + 4*7.99 + 5*0.55) / 452.70  =  1.6607

**The average is over covered sky only.** A ring of discs can enclose a hole
no telescope sees; ``hyper_fov`` drops those patches, so they don't pull the
average down. It answers *given a point the array sees at all, how many
telescopes see it*, which is why the mean and the hyper FoV are independent
numbers, not two readings of the same thing.

**The mean ignores** ``m_cut``. ``hyper_fov`` applies the cut only to the
area it returns; the patch list always includes every patch with
:math:`m \geq 1`. Passing patches from a cut run changes nothing:

.. code-block:: python

    >>> for cut in (1, 2, 3):
    ...     _, patches = array.hyper_fov(m_cut=cut)
    ...     print(cut, len(patches), round(array.multiplicity_moments(patches=patches)[0], 4))
    1 105 1.6607
    2 105 1.6607
    3 105 1.6607

This matters when reading a
:func:`~divtel.visualization.multiplicity_plot` title, which quotes the
``m_cut`` area and the mean side by side: the area respects the cut and the
mean does not.

Each telescope's disc contributes its area to every patch it covers, so
:math:`\sum_i m_i A_i` is the array's total camera solid angle, and

.. math::

   \langle m \rangle = \frac{\sum_\mathrm{telescopes} \Omega_\mathrm{camera}}
                              {\Omega_\mathrm{covered}}

the mean is that total divided by the sky actually covered. Parallel
pointing puts every camera on the same disc, so the mean is the telescope
count, with zero variance. Fully diverged, the discs are disjoint, the union
equals the sum, and the mean is one.

The mean summarizes how divergent a configuration is; it moves opposite the
area. On this array, pointed at 70 degrees:

=========  =================  ===================
``div``    hyper FoV [deg²]   mean multiplicity
=========  =================  ===================
0          46.3               16.3
0.01       100.0              7.5
0.02       169.8              4.4
0.05       452.7              1.7
=========  =================  ===================

Ten times the sky, a tenth of the telescopes on any given part of it. Which
end of that trade to pick depends on the science case.

:func:`~divtel.visualization.multiplicity_plot` draws the profile as a bar
chart, coloured to match the sky map so the two read together:

.. code-block:: python

    from divtel.visualization import display_hyper_fov, multiplicity_plot

    fig, (sky, bars) = plt.subplots(1, 2, figsize=(12, 5))
    display_hyper_fov(array, ax=sky)
    multiplicity_plot(array, ax=bars)

Both take ``m_cut``, which sets the multiplicity counted towards the area in
the title; patches and bars below it are still drawn, faded.

If you already have the patches from a
:meth:`~divtel.telescope.Array.hyper_fov` call, pass them in as ``patches=`` to
save computing the geometry twice.


Plotting
========

Two views:

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
to have that title report the stereoscopic area instead. The patches below
the cut are still drawn, faded, so you still see the full shape of the
coverage.

Both take an existing ``ax`` and return it, so they compose with whatever else
you are plotting.


Sub-arrays
==========

A real array is usually several instruments sharing a site. The La Palma
layout has four LSTs inside fifteen MSTs, each with its own camera, field of
view, and barycenter. :meth:`~divtel.telescope.Array.group_by` splits them
apart:

.. code-block:: python

    >>> array = load_array(files("divtel") / "data" / "la_palma_4LST_15MST.ecsv")
    >>> array.divergent_pointing(0.02, 70 * u.deg, 180 * u.deg)
    >>> groups = array.group_by({"LST": range(1, 5), "MST": range(5, 20)})
    >>> groups["LST"].barycenter
    <Quantity [ 0.895 , 44.8475, 44.925 ] m>
    >>> groups["LST"].hyper_fov()[0]
    <Quantity 26.83755561 deg2>
    >>> groups["MST"].hyper_fov()[0]
    <Quantity 169.76149694 deg2>

The MST figure is the whole array's coverage: the MSTs' cameras are wide enough
that they already see everything the LSTs do, so adding the LSTs back adds
multiplicity rather than area.

Each group is a full :class:`~divtel.telescope.Array`, so everything an array
can do a group can do too. The ids are
:attr:`~divtel.telescope.Telescope.id` values, not positions, which is why
layouts carry them.

Groups share the same :class:`~divtel.telescope.Telescope` objects as the
array they came from, not copies, so re-pointing the array re-points every
group too. They inherit the array's current pointing, so
:meth:`~divtel.telescope.Array.hyper_fov` works on a group right away.
Telescopes can be left out of every group; nothing requires full coverage.

To skip writing ids out, group by camera radius instead, since telescopes of
a type share one:

.. code-block:: python

    >>> {name: len(group.telescopes) for name, group in array.group_by("camera_radius").items()}
    {'1.074 m': 15, '1.05 m': 4}

:func:`~divtel.visualization.display_groups` draws them, one colour per group,
each with its own barycenter and mean pointing arrow:

.. code-block:: python

    from divtel.visualization import display_groups

    display_groups(groups)

Turn the divergence up and the two barycenter arrows swing apart: the LSTs,
sitting near the centre of the array, diverge far less than the MSTs around
them.


Pointing at a source, at a time
===============================

Everything so far works in the ground frame: an array pointed at alt 70, az 180
stays pointed there whatever the hour. Real sources rise and set, so the
direction to point at depends on where you are and when.

:class:`~divtel.observation.Observation` is a place and a time. It converts a
sky position into the alt/az pair
:meth:`~divtel.telescope.Array.divergent_pointing` takes:

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from divtel.observation import Observation
    >>> obs = Observation(site="north", time="2026-03-01T23:00:00")
    >>> crab = SkyCoord(ra=83.633 * u.deg, dec=22.015 * u.deg)
    >>> alt, az = obs.altaz_of(crab)
    >>> f"{alt:.2f}, {az:.2f}"
    '50.93 deg, 270.16 deg'
    >>> array.divergent_pointing(0.02, alt, az)

``site`` takes ``"north"`` or ``"south"`` for the two CTAO sites, built in so
this needs no network. Any other name is looked up in astropy's site
registry instead, which does need one. Pass an
:class:`~astropy.coordinates.EarthLocation` directly when the exact position
of a particular telescope matters.

An ``Observation`` is immutable. :meth:`~divtel.observation.Observation.at`
and :meth:`~divtel.observation.Observation.after` return a new one, so
watching a source across a night is a plain loop:

.. code-block:: python

    for hours in range(0, 7, 2):
        moment = obs.after(hours * u.hour)
        alt, az = moment.altaz_of(crab)
        array.divergent_pointing(0.02, alt, az)
        print(hours, array.hyper_fov()[0], array.multiplicity_moments()[0])

:attr:`~divtel.observation.Observation.sun` and
:attr:`~divtel.observation.Observation.moon` give those bodies in the same
frame, to check a time is actually observable:

.. code-block:: python

    >>> f"{obs.sun.alt:.1f}"
    '-49.9 deg'

Going the other way, :func:`~divtel.observation.pointing_coord` says where each
telescope is looking on the sky. Under divergent pointing they all look
somewhere different, so this is one coordinate per telescope:

.. code-block:: python

    >>> from divtel.observation import pointing_coord
    >>> pointing_coord(array, obs).separation(crab).to(u.deg).max()
    <Angle 4.3121837 deg>

.. note::

   divtel's azimuth is astropy's: measured from north through east, no
   conversion needed. Point an array at a source, read the telescopes back
   with :func:`~divtel.observation.pointing_coord`, and you recover the
   original position to well under an arcsecond. The test suite checks
   exactly this.


Exporting to sim_telarray
=========================

Once a pointing looks right, :meth:`~divtel.telescope.Array.export_cfg` writes
it out as a ``sim_telarray`` configuration file, so the layout can be handed to
a real shower simulation:

.. code-block:: python

    >>> array.export_cfg(outdir="configs/", tel_configs=["lst.cfg"] * 4)
    PosixPath('configs/CTA-ULTRA6-LaPalma-div0_02-az180_0-alt70_0.cfg')

The default filename records the divergence and mean pointing, so a scan
over ``div`` lands in distinct files without naming them yourself.
``tel_configs`` gives the per-telescope config each block ``#include``\ s, in
array order. Left out, it defaults to four LSTs followed by MSTs with
NectarCam, the La Palma layout the writer was built for, which is only
correct if that's your array.

``sim_telarray``'s angle conventions differ from divtel's: it takes a zenith
angle instead of altitude, and its azimuth runs the other way.
``export_cfg`` converts both on the way out.


Where to go next
================

* :doc:`examples`: the same model with sliders, running in your browser.
* :doc:`docstring`: the full API.
* `Source and issue tracker <https://github.com/cta-observatory/divtel>`_.

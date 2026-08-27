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


How well the sky is covered
===========================

The hyper FoV says how much sky the array covers.
:meth:`~divtel.telescope.Array.multiplicity_profile` says how well it covers
it -- how much of that sky is seen by exactly one telescope, by exactly two,
and so on:

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
each lying wholly inside or wholly outside every telescope's camera disc, so
each has a well-defined multiplicity :math:`m_i` and an area :math:`A_i`. The
mean and variance are those, **weighted by solid angle**:

.. math::

   \langle m \rangle = \frac{\sum_i m_i A_i}{\sum_i A_i}
   \qquad
   \mathrm{Var}(m) = \frac{\sum_i \left(m_i - \langle m \rangle\right)^2 A_i}{\sum_i A_i}

Weighting by area rather than counting patches matters: the patches are the
regions cut out by overlapping discs, so a two-telescope array makes as many
slivers as it does large regions, and an unweighted mean would let a sliver
count for as much sky as the whole field. The areas are true solid angles --
``hyper_fov`` works in an equal-area projection -- and the patches partition
the covered area exactly, so nothing is counted twice.

Working the example above by hand:

.. code-block:: text

    m :  [1       2       3      4     5   ]
    A :  [242.02  131.34  70.81  7.99  0.55]  deg²    sum = 452.70

    (1*242.02 + 2*131.34 + 3*70.81 + 4*7.99 + 5*0.55) / 452.70  =  1.6607

Two things follow from the definition that are easy to trip over.

**The average is over the covered sky, not the whole sky.** A ring of discs
encloses a hole that no telescope sees; ``hyper_fov`` drops those patches, so
they never reach this average. Read it as *given a point the array sees at all,
how many telescopes see it*. Sky the array misses does not drag the mean down,
which is exactly why the mean and the hyper FoV are two independent numbers
rather than one restating the other.

**The mean ignores** ``m_cut``. ``hyper_fov`` applies ``m_cut`` only to the
area it returns; the patch list is always every patch with :math:`m \geq 1`.
So passing patches from a cut run changes nothing:

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

There is a shortcut worth knowing. Each telescope's disc contributes its own
area to every patch it covers, so :math:`\sum_i m_i A_i` is just the total
camera solid angle of the array, and

.. math::

   \langle m \rangle = \frac{\sum_\mathrm{telescopes} \Omega_\mathrm{camera}}
                              {\Omega_\mathrm{covered}}

the total camera solid angle divided by the sky actually covered. That makes
the two limits obvious: pointed in parallel every camera lands on the same
disc, so the sum is :math:`n` times the union and the mean is the number of
telescopes with zero variance; fully diverged the discs are disjoint, the union
*equals* the sum, and the mean is one.

The mean is the summary of how divergent a configuration is, and it moves
opposite the area. On this array, pointed at 70 degrees:

=========  =================  ===================
``div``    hyper FoV [deg²]   mean multiplicity
=========  =================  ===================
0          46.3               16.3
0.01       100.0              7.5
0.02       169.8              4.4
0.05       452.7              1.7
=========  =================  ===================

Ten times the sky for a tenth of the telescopes on any part of it. Which end of
that trade you want is the question divergent pointing exists to ask.

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


Sub-arrays
==========

A real array is not one instrument but several sharing a site. On the La Palma
layout four LSTs sit inside fifteen MSTs, with different cameras, different
fields of view, and their own barycenters.
:meth:`~divtel.telescope.Array.group_by` splits the array up so each can be
looked at on its own:

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
:attr:`~divtel.telescope.Telescope.id` values, not positions -- which is why
layouts carry them.

Groups hold the same :class:`~divtel.telescope.Telescope` objects as the array
they came from rather than copies, so re-pointing the array re-points every
group with it. They also inherit the array's current pointing, so
:meth:`~divtel.telescope.Array.hyper_fov` works on a group straight away.
Telescopes may be left out of every group; nothing requires the groups to cover
the array.

If you would rather not write the ids out, grouping by camera picks the types
out on its own, since telescopes of a type share one:

.. code-block:: python

    >>> {name: len(group.telescopes) for name, group in array.group_by("camera_radius").items()}
    {'1.074 m': 15, '1.05 m': 4}

:func:`~divtel.visualization.display_groups` draws them, one colour per group,
each with its own barycenter and mean pointing arrow:

.. code-block:: python

    from divtel.visualization import display_groups

    display_groups(groups)

Turn the divergence up and the two barycenter arrows swing apart -- the LSTs,
sitting near the centre of the array, diverge far less than the MSTs around
them.


Pointing at a source, at a time
===============================

Everything so far works in the ground frame: an array pointed at alt 70, az 180
stays pointed there whatever the hour. Real sources rise and set, so the
direction to point at depends on where you are and when.

:class:`~divtel.observation.Observation` supplies that. It is a place and a
time, and it converts a sky position into the alt/az pair
:meth:`~divtel.telescope.Array.divergent_pointing` already takes:

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from divtel.observation import Observation
    >>> obs = Observation(site="north", time="2026-03-01T23:00:00")
    >>> crab = SkyCoord(ra=83.633 * u.deg, dec=22.015 * u.deg)
    >>> alt, az = obs.altaz_of(crab)
    >>> f"{alt:.2f}, {az:.2f}"
    '50.93 deg, 270.16 deg'
    >>> array.divergent_pointing(0.02, alt, az)

``site`` takes ``"north"`` or ``"south"`` for the two CTAO sites, whose
coordinates are built in so this needs no network. Any other name is looked up
in astropy's site registry, which does; an
:class:`~astropy.coordinates.EarthLocation` is taken as it is, which is what to
use when the exact position of a particular telescope matters.

An observation never changes. :meth:`~divtel.observation.Observation.at` and
:meth:`~divtel.observation.Observation.after` hand back a new one, so an array
pointed from one cannot quietly fall out of step with it. That makes watching a
source across a night a plain loop:

.. code-block:: python

    for hours in range(0, 7, 2):
        moment = obs.after(hours * u.hour)
        alt, az = moment.altaz_of(crab)
        array.divergent_pointing(0.02, alt, az)
        print(hours, array.hyper_fov()[0], array.multiplicity_moments()[0])

:attr:`~divtel.observation.Observation.sun` and
:attr:`~divtel.observation.Observation.moon` give those bodies in the same
frame, which is how you check a time is actually observable:

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

   divtel's azimuth is astropy's -- measured from north through east -- so
   nothing is converted anywhere in between. Point an array at a source and
   read the telescopes back with
   :func:`~divtel.observation.pointing_coord` and you recover the position you
   started from to well under an arcsecond, which is what the test suite
   checks.


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

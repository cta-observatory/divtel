======
divtel
======

divtel makes toy simulations of the **divergent pointing** mode for arrays of
Imaging Atmospheric Cherenkov Telescopes.

Point an array's telescopes slightly away from one another and it sees a wider
patch of sky, but fewer telescopes see any given part of it, and a shower
needs two of them to be reconstructed stereoscopically. divtel lets you set up
that trade and measure both sides of it.

.. code-block:: python

    import astropy.units as u
    from divtel.telescope import Telescope, Array

    array = Array([
        Telescope(x * u.m, y * u.m, 0 * u.m, 20 * u.m, 1 * u.m)
        for x, y in [(100, 0), (0, 100), (-100, 0), (0, -100)]
    ])

    array.divergent_pointing(0.02, 70 * u.deg, 180 * u.deg)

    covered, _ = array.hyper_fov()           # 45.96 deg2 seen at all
    stereo, _ = array.hyper_fov(m_cut=2)     # 30.27 deg2 seen by 2 or more

Start with the :doc:`guide`, or go straight to the
:doc:`interactive demo <examples>` and move the sliders.

.. toctree::
   :maxdepth: 2
   :caption: Documentation

   guide
   examples
   docstring

.. toctree::
   :maxdepth: 1
   :caption: Project

   README


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

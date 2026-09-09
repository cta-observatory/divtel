=================
API documentation
=================

Generated from the docstrings. For an introduction to what these do and how
they fit together, start with the :doc:`guide`.

``divtel.telescope``
====================

The two classes you build a simulation from.
:class:`~divtel.telescope.Telescope` is one telescope on the ground.
:class:`~divtel.telescope.Array` is a list of them, plus everything about the
array as a whole: where it points, how much sky it sees, how to write that
out.

.. automodule:: divtel.telescope
   :members:
   :undoc-members:
   :show-inheritance:

``divtel.layout``
=================

Reading an array layout from file. Layouts are ECSV tables that carry their own
units, so a camera radius may be given as a length or as an angle on the sky
without ambiguity.

.. automodule:: divtel.layout
   :members:
   :undoc-members:
   :show-inheritance:

``divtel.region``
=================

A target that is not a point: a set of sky directions with a weight on each.
Everything in :mod:`divtel.strategy` works on one, and building one from a
published sky map is the only step that needs anything beyond divtel itself.

.. automodule:: divtel.region
   :members:
   :undoc-members:
   :show-inheritance:

``divtel.strategy``
===================

Covering a weighted region: sequential tiling, simultaneous sub-arrays, and
shaped pointing that makes the array's depth follow the weight. :doc:`gw170817`
runs all of them against a real gravitational-wave localization.

.. automodule:: divtel.strategy
   :members:
   :undoc-members:
   :show-inheritance:

``divtel.pointing``
===================

The pointing geometry underneath
:meth:`~divtel.telescope.Array.divergent_pointing`. Most users won't need to
call these directly. They define the coordinate frame and the point **G**
that divergence is built around; see :ref:`the-div-parameter`.

.. automodule:: divtel.pointing
   :members:
   :undoc-members:
   :show-inheritance:

``divtel.observation``
======================

Where and when the array is observing from: the optional half that connects
the ground-frame geometry to real sources at real times. Nothing in
:mod:`divtel.telescope` depends on it.

.. automodule:: divtel.observation
   :members:
   :undoc-members:
   :show-inheritance:

``divtel.skymap``
=================

Reading a published HEALPix localization and cutting credible regions out of
it. The one part of divtel with a dependency of its own, so it lives behind an
extra::

    pip install divtel[skymap]

Precomputed GW170817 regions ship in ``divtel/data/gw170817``, so a study that
only wants those needs neither this module nor the extra.

.. automodule:: divtel.skymap
   :members:
   :undoc-members:
   :show-inheritance:

``divtel.visualization``
========================

Plotting helpers. Each takes an optional ``ax`` and returns it, so they compose
with whatever else you are drawing.

.. automodule:: divtel.visualization
   :members:
   :undoc-members:
   :show-inheritance:

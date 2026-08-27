=================
API documentation
=================

Generated from the docstrings. For an introduction to what these do and how
they fit together, start with the :doc:`guide`.

``divtel.telescope``
====================

The two classes you build a simulation from:
:class:`~divtel.telescope.Telescope` is one telescope on the ground,
:class:`~divtel.telescope.Array` is a list of them plus everything that is a
property of the array as a whole -- where it points, how much sky it sees, and
how to write that out.

.. automodule:: divtel.telescope
   :members:
   :undoc-members:
   :show-inheritance:

``divtel.pointing``
===================

The pointing geometry underneath :meth:`~divtel.telescope.Array.divergent_pointing`.
Most users will not need to call these directly, but they are what defines the
coordinate frame and the point **G** that divergence is built around -- see
:ref:`the-div-parameter`.

.. automodule:: divtel.pointing
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

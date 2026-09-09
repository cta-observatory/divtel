:hide-toc:

Interactive Examples
====================


.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - `Divergent pointing <marimo/interactive_display/index.html>`_
     - Introduction to the divergent pointing mode, with a single array. Drag the
       divergence parameter and see how the hyper field of view and mean
       multiplicity trade off.
   * - `CTAO sites <marimo/two_sites/index.html>`_
     - CTAO is two arrays: North at La Palma, 4 LSTs and 9 MSTs, and South
       at Paranal, 14 MSTs and 37 SSTs with no LSTs at all. Loads both, and
       compares what their size and camera mix buy in hyper FoV and
       multiplicity.
   * - `Observing a real source <marimo/observing_a_source/index.html>`_
     - Pick a source and a night. When is it up, when is the sun down, and what
       happens to the array's coverage as it tracks a source across the sky.
       Uses :class:`~divtel.observation.Observation`.
   * - `Choosing div <marimo/choosing_div/index.html>`_
     - The question the guide leaves open. Decide what field of view you want
       and what multiplicity you can live with, and those two fix ``div`` --
       which turns out to move as the source does.

Each is a marimo notebook under ``examples/marimo/``, editable locally with
``marimo edit examples/marimo/<name>.py``.

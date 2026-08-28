:hide-toc:

Examples
========

Divergent pointing, live
------------------------

Drag the divergence, altitude and azimuth sliders and watch the array
repoint. This runs on a Python interpreter compiled to WebAssembly: nothing
is installed, no code leaves your machine. The first load fetches the
interpreter and takes a few seconds.

.. raw:: html

    <div class="divtel-demo">
      <iframe src="marimo/interactive_display/index.html"
              title="Divergent pointing, interactive demo"
              loading="lazy"></iframe>
      <p class="divtel-demo__fallback">
        <a href="marimo/interactive_display/index.html" target="_blank" rel="noopener">Open the
        demo in its own tab</a> for the full-screen version.
      </p>
    </div>

What the panels show
--------------------

The **ground** view plots each telescope and the direction it points. The
**hyper field of view** is the sky the array actually sees: every telescope
covers a disc, and the shading counts how many telescopes see each patch.
Two or more can reconstruct a shower stereoscopically, one cannot, so
divergence buys width at the cost of depth. The area quoted in the title
counts only the part still seen by at least two.

The cameras here are about 5.7 degrees across, while divergence swings the
pointings by tens of degrees, so the discs come apart within the first few
hundredths of the slider. Push it past roughly 0.08 and the array covers
four separate patches of sky with no stereoscopic overlap at all.

The sky panel is drawn in an equal-area projection centred on the array's
mean pointing: its axes are degrees of offset from that direction, and the
quoted area is a true solid angle wherever the array points.

The demo is built from ``examples/marimo/interactive_display.py``. To run or
edit it locally::

    marimo edit examples/marimo/interactive_display.py

A Jupyter version lives in
``examples/notebooks/interactive_display.ipynb``, but isn't part of this
site: its ``ipywidgets`` sliders need a Python kernel, which a static page
doesn't have. Run it in a local Jupyter session instead.


Tutorials
---------

Four worked notebooks, running in your browser the same way. They pick up
where the :doc:`guide` leaves off: the guide says what the API does, these
ask what to do with it.

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - `Two sites <marimo/two_sites/index.html>`_
     - CTAO is two arrays: North at La Palma, 4 LSTs and 9 MSTs, and South
       at Paranal, 14 MSTs and 37 SSTs with no LSTs at all. Loads both, and
       compares what their size and camera mix buy in hyper FoV and
       multiplicity.
   * - `One site, two instruments <marimo/sub_arrays/index.html>`_
     - The CTAO North layout is four LSTs inside nine MSTs, with different
       cameras. Why the hyper FoV of the whole array is the MSTs' number, why
       one ``div`` means two different angles, and what
       :meth:`~divtel.telescope.Array.group_by` is for.
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

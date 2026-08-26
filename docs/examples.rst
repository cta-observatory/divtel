:hide-toc:

Examples
========

Divergent pointing, live
------------------------

Drag the divergence, altitude and azimuth sliders and watch the array
repoint. This runs on a Python interpreter compiled to WebAssembly, so it
executes in your browser: nothing is installed and no code leaves your
machine. The first load fetches the interpreter and takes a few seconds.

.. raw:: html

    <div class="divtel-demo">
      <iframe src="marimo/index.html"
              title="Divergent pointing, interactive demo"
              loading="lazy"></iframe>
      <p class="divtel-demo__fallback">
        <a href="marimo/index.html" target="_blank" rel="noopener">Open the
        demo in its own tab</a> for the full-screen version.
      </p>
    </div>

What the panels show
--------------------

The **ground** view plots each telescope and the direction it points. The
**hyper field of view** is the sky the array actually sees: every telescope
covers a disc, and the shading counts how many telescopes see each patch.
Two or more can reconstruct a shower stereoscopically, one cannot -- so
divergence buys width at the cost of depth, and the area quoted in the title
counts only the part still seen by at least two.

The cameras in this example are about 5.7 degrees across while divergence
swings the pointings by tens of degrees, so the discs come apart within the
first few hundredths of the slider. Push it past roughly 0.08 and the array
covers four separate patches of sky with no stereoscopic overlap at all.

The sky panel is drawn in an equal-area projection centred on the array's
mean pointing, so its axes are degrees of offset from that direction and the
quoted area is a true solid angle, whatever the array is pointing at.

The demo is built from ``examples/marimo/interactive_display.py``. To run or
edit it locally::

    marimo edit examples/marimo/interactive_display.py

A Jupyter version of the same example lives in
``examples/notebooks/interactive_display.ipynb``. It is not part of this site:
its ``ipywidgets`` sliders need a Python kernel, which a static page has not
got, so it would render here as controls that cannot move. Run it in a local
Jupyter session instead.

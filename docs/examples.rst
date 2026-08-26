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

The demo is built from ``examples/marimo/interactive_display.py``. To run or
edit it locally::

    marimo edit examples/marimo/interactive_display.py

A Jupyter version of the same example lives in
``examples/notebooks/interactive_display.ipynb``. It is not part of this site:
its ``ipywidgets`` sliders need a Python kernel, which a static page has not
got, so it would render here as controls that cannot move. Run it in a local
Jupyter session instead.

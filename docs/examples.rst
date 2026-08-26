Examples
========

Divergent pointing, live
------------------------

Drag the divergence, altitude and azimuth sliders and watch the array
repoint. This runs on a Python interpreter compiled to WebAssembly, so it
executes in your browser: nothing is installed and no code leaves your
machine. The first load fetches the interpreter and takes a few seconds.

.. raw:: html

    <style>
      /* The Read the Docs content column is far narrower than the demo wants,
         so let this one block escape it on screens wide enough to spare the
         room. */
      .divtel-demo {
        margin: 1.5em 0;
      }
      .divtel-demo iframe {
        width: 100%;
        height: 1250px;
        border: 1px solid #e1e4e5;
        border-radius: 4px;
        background: #fff;
      }
      @media (min-width: 1400px) {
        .divtel-demo {
          width: 1150px;
          margin-left: calc((100% - 1150px) / 2);
        }
      }
      .divtel-demo__fallback {
        font-size: 90%;
        margin-top: 0.5em;
      }
    </style>

    <div class="divtel-demo">
      <iframe src="marimo/index.html"
              title="Divergent pointing, interactive demo"
              loading="lazy"></iframe>
      <p class="divtel-demo__fallback">
        Cramped in here?
        <a href="marimo/index.html" target="_blank" rel="noopener">Open the
        demo in its own tab</a>.
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

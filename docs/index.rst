:hide-toc:

======
divtel
======

divtel makes toy simulations of the **divergent pointing** mode for arrays of
Imaging Atmospheric Cherenkov Telescopes.

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


The **ground** view plots each telescope and the direction it points. The
**hyper field of view** is the sky the array actually sees: every telescope
covers a disc, and the shading counts how many telescopes see each patch.
Two or more can reconstruct a shower stereoscopically, one cannot, so
divergence buys width at the cost of depth. The area quoted in the title
counts only the part still seen by at least two.


.. toctree::
   :maxdepth: 2
   :caption: Documentation

   guide
   examples
   docstring

.. toctree::
   :maxdepth: 1
   :caption: Studies

   introduction
   definitions
   ceiling
   tracking

.. toctree::
   :maxdepth: 1
   :caption: Project

   README


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

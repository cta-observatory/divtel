============
Introduction
============

An imaging atmospheric Cherenkov telescope sees a few degrees of sky. Point
every telescope in an array the same way and their fields of view land on
top of one another, so the array sees what its widest camera sees and no
more: CTAO-North's widest camera has an angular radius of 3.8°, giving the
array 46 deg²; CTAO-South's is 4.4°, giving 61 deg². Fifty-one telescopes do
not watch fifty-one times more sky than one of them — pointed conventionally,
adding telescopes buys sensitivity, not sky.

Several observing programmes CTAO plans need more sky than that at one
instant. The Extragalactic Survey proposes to cover a quarter of the sky in
a fixed time budget set by how much sky one pointing covers. A
gravitational-wave alert arrives with a 90% credible region anywhere from
ten to a thousand square degrees, and these regions are rarely compact: the
first sky map circulated for GW170817 covered 187 deg², but its containment
radius, the angular radius of the smallest cone enclosing it, was 62.1° —
against a camera radius of about 4°, the shortfall is in reach, not in
collecting area.

There are three ways to close that gap. **Tiling** keeps the telescopes
together and visits the region in pieces, spending time instead of sky.
**Divergent pointing** tilts the telescopes apart so their fields of view
overlap only partially, widening the instantaneous footprint at the cost of
the number of telescopes seeing any given direction. **Splitting the array**
points independent groups at different parts of the region at once. All
three are geometric decisions, and all three are constrained by the same
fact: a shower must be recorded by at least two telescopes to be
reconstructed stereoscopically, so sky seen once is not sky observed.

A large part of that trade-off does not need a shower simulation to
evaluate. Given a layout and a set of pointings, whether a sky direction
falls inside a camera is a geometric fact. The number of telescopes
containing a given direction — the **multiplicity** — is therefore exact,
and so is everything built from it: the combined field of view, the part
of it that is stereoscopic, and the divergence at which stereoscopic
coverage stops improving. What geometry cannot supply is the map from
multiplicity to sensitivity, which is a separate measurement, not an
assumption this package makes.

``divtel`` is the open-source package that computes these geometric
quantities for an arbitrary array under arbitrary per-telescope pointing.
This section of the documentation follows the structure of the paper that
introduces it, T. Vuillaume, *Planning divergent and non-traditional
pointing for imaging atmospheric Cherenkov telescope arrays with divtel*,
as far as the divergence ceiling:

* :doc:`definitions` sets out the quantities ``divtel`` computes — the
  hyper field of view, the multiplicity, and the mean multiplicity — since
  everything that follows is written in them, and shows what an array and
  its coverage look like.
* :doc:`capabilities` points at a real source and tracks it across a
  night, interactively.
* :doc:`ceiling` answers the question the rest of this section builds
  towards: how far an array can usefully be spread, and what that is worth
  for CTAO-North and CTAO-South.

The paper's later sections — tiling against divergence, shaped pointing
that follows a sky map, and matching a configuration to a science case —
are not covered here yet.

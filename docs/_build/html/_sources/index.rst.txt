.. raw:: html

  <p>
    <a href="https://badge.fury.io/py/polycleaver"><img alt="pypi badge" src="https://badge.fury.io/py/polycleaver.svg"/></a>
  </p>


PolyCleaver · an ionic surface generation package
=================================================

:program:`PolyCleaver` PolyCleaver is a Python-based package that generates 
high quality vacuum-containing surfaces from bulk structures of mineral structures
characterised as ionic compounds with polyatomic anions (e.g., Mg\ :sub:`2`\SiO\ :sub:`4`).
These include silicates, sulfides, carbonates, sulfates and phosphates, among others.

.. card:: Access the web app
   :link: https://polycleaver.streamlit.app
   :class-card: sd-text-white sd-text-center sd-bg-info

This algorithm is built around the pymatgen library, allowing for a high degree of
customization and future enhancement for other ionic compounds.
Surfaces generated using the PolyCleaver algorithm are:

- Non-polar, allowing accurate surface reactivity calculations.
- Stoichiometric with respect to their bulk composition, maintaining per-atom charges.
- Due to the high energy nature of the bonds in the covalent units forming the polyatomic
  anions (e.g. SiO\ :sub:`4`\ :sup:`2-`\), cleave is carried out maintaining all covalent bonds.

This algorithm detects all structural parameters of the bulk automatically
(e.g. identification of species, clustering of covalent units, calculation of
coordination numbers) and performs a sub-set of cuts. These slabs are then corrected
using a series of symmetry and geometry operations to generate the final structures.
Geometrical parameters of the slabs (e.g. thickness, number of undercoordinated cations
on the topmost layers) are easily accessible, facilitating an unsupervised high-throughput
generation of surface slabs with any given set of Miller indices.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

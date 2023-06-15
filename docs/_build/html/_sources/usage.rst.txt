.. _usage:
.. index:: Usage

Usage
*****

This package has been implemented in a Streamlit-based web app, accessible :bdg-link-primary:`here <https://polycleaver.streamlit.app>`.

Ionic slabs with polyatomic anions
==================================

Given a bulk structure of an ionic compound with a polyatomic anionic subunit
(e.g., Mg\ :sub:`2`\SiO\ :sub:`4`), :program:`PolyCleaver` is able to generate
all surfaces from a list of Miller indices as follows:

.. code-block:: python

    from polycleaver.core import mnx
    miller_indices_list = [(0,0,1), (1,0,1)]
    slabs = mnx.generate_mnx_slabs('bulk_mg2sio4.cif', miller_indices_list, save=True)

This will generate all non-polar terminations of the selected bulk with the specified
Miller indices. The approximate thickness and vacuum of the slab (by default, 15 Å each)
can be specified as follows:

.. code-block:: python
   
    slabs = mnx.generate_mnx_slabs('bulk_mg2sio4.cif', miller_indices_list, thickness=15, vacuum=15, save=True)

This will save all structures in cif files, with an integer index plus the Miller three
indices of the generated plane (e.g., 1_001.cif). Alternatively, if planes with all 
`hkl` indices up to a maximum Miller index of `n` are to be modelled, instead of a Miller
indices list, the user can specify a single integer corresponding to `n`, and the package will
construct a list containing all non-equivalent Miller indices based on the bulk symmetry:

.. code-block:: python

    from polycleaver.core import mnx
    max_index = 2
    slabs = mnx.generate_mnx_slabs('bulk_mg2sio4.cif', max_index, save=True)

Visualization of slabs
----------------------

The slabs can be visualized using ASE as follows:

.. code-block:: python
   
    from pymatgen.io.ase import AseAtomsAdaptor
    from ase.visualize import view

    view([AseAtomsAdaptor.get_atoms(slab.atoms) for slab in slabs])

Calculation of slab parameters
------------------------------

After generation of surfaces, :program:`PolyCleaver` can be used to display relevant surface
parameters (only working for ionic surfaces with polyatomic anions, generated with the
`generate_mnx_slabs` module):

.. code-block:: python
   
    for slab in slabs:
        print(slab.thickness)    # Prints the final thickness of the reconstructed slabs
        print(slab.cations)      # Prints list of all cationic atoms (e.g., Mg2+ in Mg2SiO4).
        print(slab.anions)       # Prints list of all anionic atoms (e.g., O2- in Mg2SiO4).
        print(slab.centers)      # Prints list of the center atoms of all polyatomic anions (e.g., Si4+ in Mg2SiO4).
        print(slab.undercoordinated_sites(slab.cations)) # Prints all undercoordinated cations at both sides of the slab.

This allows for rapid selection of most interesting slabs for reactivity or adsorption studies.

Ionic slabs with monoatomic anions
==================================

Alternatively, this package can also generate bare slabs from ionic compounds with
monoatomic ions such as FeS\ :sub:`2` in a similar fashion, using the `polycleaver.core.mx` module:

.. code-block:: python
   
    from polycleaver.core import mx
    miller_indices_list = [(0,0,1), (1,0,1)]
    slabs = mx.generate_mx_slabs('bulk_fes2.cif', miller_indices_list, save=True)


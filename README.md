# PolyCleaver

A repository for the generation of non-polar, neutral surfaces from ionic compounds with polyatomic anions.

PolyCleaver is a robust Python-based algorithm that allows for the generation of high quality vacuum-containing surfaces from bulk structures of mineral structures characterised as **ionic compounds with polyatomic anions** (e.g., Mg<sub>2</sub>SiO<sub>4</sub>). These include **silicates**, **sulfides**, **carbonates**, **oxides**, **halides**, **sulfates** and **phosphates**, among others. Similar structures, such as triatomic **perovskites**, are also supported. This algorithm is built around the *pymatgen* library, allowing for a high degree of customization and future enhancement for other ionic compounds. Surfaces generated using the PolyCleaver algorithm are:

1. **Non-polar**, allowing accurate surface reactivity calculations.
2. **Stoichiometric** with respect to their bulk composition, maintaining per-atom charges.
3. Due to the high energy nature of the bonds in the covalent units forming the polyatomic anions (e.g. SiO<sub>4</sub><sup>2-</sup>), cleave is carried out **maintaining all covalent bonds**. 

This algorithm detects all structural parameters of the bulk automatically (e.g. identification of species, clustering of covalent units, calculation of coordination numbers) and performs a sub-set of cuts using pymatgen's SlabGenerator class. These slabs are then corrected using a series of symmetry and geometry operations to generate the final structures. Geometrical parameters of the slabs (e.g. thickness, number of undercoordinated cations on the topmost layers) are easily accessible, facilitating an unsupervised high-throughput generation of surface slabs with any given set of Miller indices.


## Usage/Examples

Starting from a given bulk structure as a pymatgen.core.Structure object:

```python
# We can load a structure from a cif file as follows:

from pymatgen.core import Structure

bulk = Structure.from_file('bulk.cif')

# Or, alternatively, use the MaterialsProject API:

from pymatgen.ext.matproj import MPRester

with MPRester("API_Key (obtainable in the MaterialsProject website)") as m:
    bulk = m.get_structure_by_material_id("mp-3164")

# Generating all possible slabs is straightforward:

from PolyCleaver import *

hkl = [(0,0,1), (0,1,0)] # We can specify the set of Miller indices to generate the surfaces,
hkl = 1                  # or we can specify the maximum Miller index that we want to analyse, 
                         # and the algorithm will investigate all non-equivalent Miller indices
                         # using pymatgen's get_symmetrically_distinct_miller_indices function.
                           
slabs = generate_slabs(bulk, hkl)
```

All slabs will then be generated and stored as a list of SlabUnit objects, which are a wrap of the Slab object of pymatgen containing a set of additional structural parameters. For instance, for any slab in "slabs":

- slab.atoms: accesses the Slab object of pymatgen, which we can then save as follows:

            slab.atoms.to_file('file.cif')

- slab.anions / slab.cations / slab.centers: outputs all anions, cations and cationic centers of the polyatomic anions as a list of sites.

- slab.undercoordinated_sites(sites): gives an integer with the number of undercoordinated sites in the slab from a list of sites. This can be combined   with the previous attribute to obtain the number of undercoordinated surface cations in the final slabs:

            slab.undercoordinated_sites(slab.cations)

- slab.thickness: gives the thickness of the slab, in Å.

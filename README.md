# PolyCleaver

A repository for the generation of non-polar, neutral surfaces from ionic compounds with polyatomic anions.

PolyCleaver is a Python-based pipeline that allows for the generation of high quality vacuum-containing surfaces from bulk structures of mineral structures characterised as **ionic compounds with polyatomic anions** (e.g., Mg<sub>2</sub>SiO<sub>4</sub>). These include **silicates**, **sulfides**, **carbonates**, **sulfates** and **phosphates**. This algorithm is built around the *pymatgen* library, allowing for a high degree of customization and future enhancement for other ionic compounds. Surfaces generated using the PolyCleaver algorithm are:

1. **Non-polar**, allowing accurate surface reactivity calculations.
2. **Stoichiometric** with respect to their bulk composition, maintaining per-atom charges.
3. Due to the high energy nature of the bonds in the covalent units forming the polyatomic anions (e.g. SiO<sub>4</sub><sup>2-</sup>), cleave is carried out **maintaining all covalent bonds**. 

This algorithm detects all structural parameters of the bulk automatically (e.g. identification of species, clustering of covalent units, calculation of coordination numbers) and performs a sub-set of cuts using pymatgen's SlabGenerator class. These slabs are then corrected using a series of symmetry and geometry operations to generate the final structures. Geometrical parameters of the slabs (e.g. thickness, number of undercoordinated cations on the topmost layers) are easily accessible, facilitating an unsupervised high-throughput generation of surface slabs with any given set of Miller indices.

For Installation, Usage and Examples, see [PolyCleaver's documentation](https://polycleaver.readthedocs.io/en/latest/index.html).

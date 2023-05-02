from pymatgen.core import Structure
from ase.visualize import view
from pymatgen.io.ase import AseAtomsAdaptor as ac
from tools import remove_equivalent_slabs, get_initial_slabs, set_site_attributes
import copy
import os
os.chdir('/Users/eric/Projects/PolyCleaver/polycleaver/core')
import numpy as np

bulk = Structure.from_file('../../examples/bulk_fes2.cif')
bulk.add_oxidation_state_by_guess()
hkl = 1
thickness = 15
vacuum = 15
initial_slabs = get_initial_slabs(bulk, hkl, thickness, vacuum)

final_slabs = [slab for slab in initial_slabs if slab.is_symmetric()]
for slab in final_slabs:
    print(slab.miller_index)
final_slabs[2].is_polar()
set_site_attributes(final_slabs[2])

[[atom.element, atom.coordination_number] for atom in final_slabs[2]]

initial_slabs[0][15]


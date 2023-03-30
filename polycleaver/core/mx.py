from pymatgen.core import Structure
import tools

bulk = Structure.from_file

initial_slabs = get_initial_slabs(bulk, hkl, thickness, vacuum)
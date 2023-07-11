import numpy as np
from . import tools
import os

def generate_mx_slabs(bulk_str, hkl, thickness=15, vacuum=15, save=True):
    """
    Generates non-polar, stoichiometric slabs from a given ionic crystal

    Args:
        bulk:       pymatgen.core.structure.Structure object of the bulk
                    of which slabs are to be generated.
        hkl:        list of Miller indices as tuples (e.g., [(0,0,1), (0,1,0)])
                    of planes parallel to surfaces to generate. Alternatively,
                    maximum integer of the Miller index to analyse (e.g., if hkl=1,
                    all symmetrically distinct Miller indices up to (1,1,1) will be
                    considered, as determined by pymatgen's
                    'get_symmetrically_distinct_miller_indices' function).
        thickness:  thickness of the preliminary slabs generated. Note that
                    this thickness will certainly vary in the final, corrected
                    slabs.
        vacuum:     vacuum of the preliminary slabs in the direction normal to
                    the surface. Note that this thickness will certainly vary
                    in the final, corrected slabs.

    Returns:
        (list): SlabUnit objects of corrected slabs with all required parameters
                for characterisation.
    """
    bulk = tools.read_bulk(bulk_str)
    bulk.add_oxidation_state_by_guess()
    initial_slabs = tools.get_initial_slabs(bulk, hkl, thickness, vacuum)
    slabs = [slab for slab in initial_slabs if not slab.is_polar()]
    print(f'{len(slabs)} non-polar slabs generated.')

    if save:
        for number, slab in enumerate(slabs):
            miller = f'{slab.atoms.miller_index[0]}{slab.atoms.miller_index[1]}{slab.atoms.miller_index[2]}'
            slab.atoms.to(f'{number}-{miller}', 'cif')

    return slabs

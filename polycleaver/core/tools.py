from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.surface import SlabGenerator, get_symmetrically_distinct_miller_indices

def remove_equivalent_slabs(slablist):
    """
    Often, the sanitizing operations performed in the 'generate_slabs' function
    converge to symmetrically equivalent slabs. This function analyses a
    list containing SlabUnit objects and removes equivalent slabs.

    Args:
        slablist:   list of SlabUnit objects from which equivalent slabs
                    are to be removed.
    """
    for index, slab in enumerate(slablist):
        equivalences = [StructureMatcher().fit(slablist[index].atoms, slablist[_index].atoms)
                        for _index in range(len(slablist))
                            if _index != index]
        if any(equivalences):
            _equivalences = [_index for _index in range(len(equivalences))
                                if equivalences[_index]
                                and _index != index
                                and slablist[index].atoms.miller_index == slablist[_index].atoms.miller_index]
            for _index in sorted(_equivalences, reverse=True):
                del slablist[_index]

def get_initial_slabs(bulk, hkl_list, thickness, vacuum, tolerance=.01):
    """
    Desc.
    """
    all_slabs = []
    print('Generating preliminary slabs...')
    if isinstance(hkl, int):
        hkl_list = get_symmetrically_distinct_miller_indices(bulk, hkl)
    if isinstance(hkl, list):
        hkl_list = hkl
    print(f'Miller indices that will be analysed: {", ".join([str(hkl) for hkl in hkl_list])}')
    for hkl in hkl_list:
        slabgen = SlabGenerator(bulk, hkl, thickness, vacuum, center_slab=True)
        all_slabs.extend(slabgen.get_slabs())

    valid_slabs = []
    for slab in all_slabs:
        if not slab.is_polar(tol_dipole_per_unit_area=tolerance):
            valid_slabs.append(slab)
        else:
            valid_slabs.append(slab)
            slab_tk_l = slab.get_tasker2_slabs(tol=tolerance)
            for slab_tk in slab_tk_l:
                valid_slabs.append(slab_tk)

    print(f'{len(valid_slabs)} preliminary slabs generated.')

    return valid_slabs
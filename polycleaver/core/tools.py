from pymatgen.analysis.structure_matcher import StructureMatcher

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
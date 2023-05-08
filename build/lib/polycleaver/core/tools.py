from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.surface import SlabGenerator, get_symmetrically_distinct_miller_indices
import numpy as np

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

def get_initial_slabs(bulk, hkl, thickness, vacuum, tolerance=.01):
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

def set_site_attributes(structure):
    """
    Sets the attributes of the sites present inside a structure.

    Args:
        structure: site pymatgen structure (Slab, Structure...) of which
                    sites are to be analysed.
    """

    def coordination_number(site, structure):
        """
        Determines the coordination number of any given periodic site
        inside a structure.

        Args:
            site: pymatgen.core.sites.PeriodicSite site to analyse.
            structure: pymatgen Structure / Slab containing the site.

        Returns:
            (int): coordination number of the selected site.
        """
        try:
            distance_nn = min(
                                [structure.get_distance(site.index, atom.index)
                                    for atom in structure.get_neighbors(site, 3.5)]
                                ) + 0.2
            neighborlist = structure.get_neighbors(site, distance_nn)
            coordination_number = len(neighborlist)
        except ValueError:
            coordination_number = 0
        return coordination_number

    def cluster(site, structure):
        """
        Performs clustering analyses of the surrounding species from a
        single site center.

        Args:
            site: pymatgen.core.sites.PeriodicSite site to analyse.
            structure: pymatgen Structure / Slab containing the site.

        Returns:
            (list): PeriodicSite objects including the center atom
            and its nearest neighbors.
        """
        cluster = []
        cluster.extend(structure.get_neighbors(site, 2.3))
        cluster.append(site)
        return np.array(cluster, dtype=object)

    #####
    ## Attributes of all pymatgen.core.sites.PeriodicSite objects
    ## conforming the structure are update here.
    ####
    for site in structure:
        site.index = structure.index(site)
        site.cluster = cluster(site, structure)
        site.coordination_number = coordination_number(site, structure)
        site.element = site.specie.element.symbol
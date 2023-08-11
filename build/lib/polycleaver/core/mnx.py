from ase.visualize import view
from pymatgen.io.ase import AseAtomsAdaptor as ac
import numpy as np
from . import tools
import copy, sys
from .miscellaneous import load_bar

class BulkUnit():
    """
    This class analyses the parameters of a monometallic
    ionic-covalent bulk structure (i.e., Mg2SiO4, CaCO3)
    based on geometrical and oxidation state data.
    """

    def __init__(self, bulk):
        """
        Makes a BulkUnit structure, wrapped around a Structure
        pymatgen class, containing geometrical parameters of
        the cations (e.g. Mg2+), anions (e.g. O2-) and
        anion centers (e.g. Si4-). A pymatgen Structure object
        can be obtained by calling the 'atoms' attribute.
        The element forming anions, centers and cations can be
        obtained with the 'anion', 'center' and 'cation' attributes,
        respectively.

        Args:
            bulk: pymatgen.core.structure.Structure object.
        """
        self.atoms = bulk
        self.bulk = self
        tools.set_site_attributes(self.atoms)
        global anion_str, center_str, cations_strs
        anion_str = self.bulk.anions[0].element
        center_str = self.bulk.centers[0].element
        cations_strs = [ atom.element for atom in list({object_.species_string: object_ for object_ in self.bulk.cations}.values())]

    @property
    def anions(self):
        """
        Define the anions in the structure, by detecting whether negative
        oxidation states are present.

        Raises:
            ValueError: if no anions are found in the structure.

        Returns:
            (list): PeriodicSite objects of all anions.
        """
        sites = np.array(self.atoms.sites)
        mask = np.vectorize(lambda site: site.specie.oxi_state < 0)(sites)
        anions = sites[mask]
        if len(anions) == 0:
            raise ValueError('No anions could be found.')
        return anions

    @property
    def centers(self):
        """
        Define the anion centers in the structure. These are detected by
        obtaining the nearest neighbours to the anion species.

        Raises:
            ValueError: if no polyatomic anions are found in the structure.

        Returns:
            (list): PeriodicSite objects of all covalent anionic unit centers.
        """
        sites = np.array(self.atoms.sites)
        species = list({object_.species_string: object_ for object_ in self.atoms}.values())
        center = sorted(species, key=lambda x: x.specie.oxi_state)[-1].specie
        if center.oxi_state < 0:
            raise ValueError('No polyatomic anions could be found.')
        mask = np.vectorize(lambda site: site.specie == center)(sites)
        centers = sites[mask]
        return centers

    @property
    def cations(self):
        """
        Define the metallic cation by determining species with positive
        oxidation states not matching centers' species.

        Raises:
            ValueError: if no cations are found in the structure.

        Returns:
            (list): PeriodicSite objects of all cations.
        """
        sites = np.array(self.atoms.sites)
        species = list({object_.species_string: object_ for object_ in self.atoms}.values())
        cation_species = []
        for specie in species:
            if specie not in self.bulk.centers and specie not in self.bulk.anions:
                cation_species.append(specie.specie)
        if len(cation_species) == 0:
            raise ValueError('No cations could be found.')
        mask = np.vectorize(lambda site: site.specie in cation_species)(sites)
        cations = sites[mask]
        return cations

class SlabUnit(BulkUnit):
    """
    Inherited BulkUnit class, adapted for the generated slabs.
    Adds additional parameters, as well as analysis functions
    to perform stoichiometry and polarity corrections.
    """

    def __init__(self, slab, bulk_obj):
        """
        Makes SlabUnit object with additional slab-oriented parameters.
        A pymatgen Slab object can be obtained by calling the 'atoms'
        attribute. The original bulk structure is stored in the
        'bulk' attribute.
        Args:
            slab: pymatgen.core.surface.Slab object.
            bulk_obj: BulkUnit object.
        """
        super().__init__(bulk_obj.atoms)
        self.atoms = slab
        self.bulk = bulk_obj
        tools.set_site_attributes(self.atoms) # Update attributes of periodic sites.

    def undercoordinated_sites(self, sites):
        """
        Retrieves number of undercoordinated sites of a surface. The coordination
        number of the saturated site is defined as the maximum coordination number
        of all the sites in the list; therefore, this function is to be used with
        all species of one type in the structure at the same time (by calling the
        'centers/anions/cations' functions in the BulkUnit and SlabUnit classes)
        Args:
            sites: list of PeriodicSite objects to investigate.
        Returns:
            (int): number of undercoordinated sites.
        """
        coord_mask = np.vectorize(lambda site: site.element == sites[0].element)(self.bulk.atoms)
        coord_no_bulk = np.vectorize(lambda site: site.coordination_number)(self.bulk.atoms)[coord_mask][0]
        undercoord_mask = np.vectorize(lambda site: site.coordination_number < coord_no_bulk)(sites)
        return sites[undercoord_mask]
    
    @property
    def thickness(self):
        """
        Determines the thickness of the slab based on the maximum and minimum
        Z cartesian positions of all atoms.
        Returns:
            (float): thickness of the slab in angstroms.
        """
        # next(atoms[-1].z - atoms[0].z for atoms in sorted(self.atoms, key=lambda x: x.z))
        atomslist = sorted(self.atoms, key=lambda x: x.z)
        return atomslist[-1].z - atomslist[0].z

    @property
    def lone_anions(self):
        """
        Analyses the structure to find the anions which are not bound todo
        their center atom (e.g., in the case of silicates, O atoms which are
        not bound to Si).
        Returns:
            (list): PeriodicSite objects of lone anions.
        """        
        def convert_to_strings(cluster):
            return np.vectorize(lambda site: site.specie.element.symbol, otypes=[object])(cluster)
        
        cluster_atoms = np.vectorize(lambda site: site.cluster)(self.anions)
        cluster_elements = np.vectorize(convert_to_strings)(cluster_atoms)
        mask_lone_anions = np.vectorize(lambda cluster: center_str not in cluster)(cluster_elements)
        return self.anions[mask_lone_anions]

    def there_are_cations(self):
        return all([cation in self.atoms.formula for cation in cations_strs])

    @property
    def clusters_to_remove(self):
        """
        Detects undercoordinated anionic centers and retrieves their nearest
        neighbors.
        Returns:
            (list): PeriodicSite objects of undercoordinated anionic centers,
            as well as their nearest anion neighbors.
        """
        try:
            atoms_to_remove = np.vectorize(lambda site: site.cluster)(self.undercoordinated_sites(self.centers))
            atoms_to_remove = np.concatenate(atoms_to_remove).ravel()
        except ValueError:
            atoms_to_remove = np.array([], dtype=object)
        return atoms_to_remove

    def remove_sites(self, sites):
        """
        Wrapping of the remove_sites function of pymatgen, this function
        removes a set of sites inside the SlabUnit class, without generating
        a copy of the structure and maintaining site attributes.
        Args:
            sites: list of PeriodicSite objects of sites which are to be
            removed from the main SlabUnit object.
        """
        if len(sites) > 0:
            template = self.atoms.copy()
            sites_to_remove = np.vectorize(lambda site:
                                                np.where(np.vectorize(lambda atom: atom.is_periodic_image(site), otypes=[object])(template))[0], 
                                        otypes=[object])(np.array(sites, dtype=object))
            sites_to_remove = np.squeeze(
                                  np.concatenate(
                                    [a for a in sites_to_remove.flatten() if len(a) > 0]
                                                )
                                        ) if len(sites_to_remove) > 1 else sites_to_remove[0]
            if sites_to_remove.ndim == 0:
                sites_to_remove = np.resize(sites_to_remove, (1, 1))
            template.remove_sites(list(sites_to_remove))
            for site in [atom for atom in self.atoms if atom not in template]:
                self.atoms.remove(site)

    def remove_element(self, elements):
        """
        Similarly to the 'remove_sites' function, this function removes
        all sites of a certain element.
        Args:
            element: element to be removed, as string (e.g., 'Fe').
        """
        mask_elements = np.vectorize(lambda site: site.element in elements)(self.atoms)
        atoms_list = np.array(self.atoms)[mask_elements]
        for atom in atoms_list:
            self.atoms.remove(atom)

    def top_site(self, element):
        """
        Finds the top site of a slab.
        Args:
            element: element of site to be investigated, as string (e.g., 'Fe').
        """
        mask_elements = np.vectorize(lambda atom: atom.element == element)(self.atoms)
        atoms_list = list(np.array(self.atoms)[mask_elements])
        atoms_list.sort(key=lambda x: x.z)
        return atoms_list[-1]

    def bottom_site(self, element):
        """
        Finds the bottom site of a slab.
        Args:
            element: element of site to be investigated, as string (e.g., 'Fe').
        """
        mask_elements = np.vectorize(lambda atom: atom.element == element)(self.atoms)
        atoms_list = list(np.array(self.atoms)[mask_elements])
        atoms_list.sort(key=lambda x: x.z)
        return atoms_list[0]

    def depolarize_anions(self):
        """
        Removes polyatomic anions from one side or the other
        from a structure without cations until non-polarity is achieved.
        In addition, returns approach value. By default, polyatomic
        anion clusters are removed from the top; if non-polarity cannot
        be achieved through this method, the structure is regenerated
        and the algorithm repeats from the bottom, changing
        'approach_value' to 'bot'.
        Returns:
            (str): side from which polyatomic anions have been removed,
            which determines the side from which cations should be
            removed in subsequent geometrical operations to
            achieve non-polarity of the whole structure.
        """
        template = copy.deepcopy(self)
        sites_to_remove = []
        approach_value = 'top'
        while (template.atoms.is_polar(tol_dipole_per_unit_area=0.01) and
               center_str in template.atoms.formula or
               len(template.centers)%2 != len(template.bulk.centers)%2 if (center_str in template.atoms.formula
                                                                       and anion_str in template.atoms.formula)
                                                                       else False):
                _sites_to_remove = template.top_site(center_str).cluster
                template.remove_sites(_sites_to_remove)
                sites_to_remove.extend(_sites_to_remove)
        if not center_str in template.atoms.formula:
            template = copy.deepcopy(self)
            sites_to_remove = []
            approach_value = 'bot'
            sites_to_remove = []
            while (template.atoms.is_polar(tol_dipole_per_unit_area=0.01) and
                   center_str in template.atoms.formula or
                   len(template.centers)%2 != len(template.bulk.centers)%2 if (center_str in template.atoms.formula
                                                                           and anion_str in template.atoms.formula)
                                                                           else False):
                _sites_to_remove = template.bottom_site(center_str).cluster
                template.remove_sites(_sites_to_remove)
                sites_to_remove.extend(_sites_to_remove)
        self.remove_sites(sites_to_remove) if (template.atoms.formula != '' and len(sites_to_remove) > 0) else None
        return approach_value

    def depolarize_cations(self, topbot=None):
        """
        Once non-polarity is achieved in the scaffold containing only
        polyatomic anions, cations are removed from one side, as
        determined by the 'topbot' parameter, until non-polarity is
        achieved (or until bulk stroichiometry is achieved in the slab).
        Args:
            topbot: side from which polyatomic anions have been removed.
        """
        atomlist = sorted(self.cations, key=lambda x: x.z)
        a = 0
        ratio_bulk = len(self.bulk.cations)/len(self.bulk.centers)
        while (self.atoms.is_polar(tol_dipole_per_unit_area=0.01) and
               len(self.cations)/len(self.centers) > ratio_bulk):
            a += -1 if topbot == 'top' else 0
            self.remove_sites([atomlist[a]])
            a += 1 if topbot == 'bot' else 0

    def stoichiometrize(self):
        """
        Removes cations from a non-polar slab until composition
        matches the bulk stoichiometry.
        """
        atomlist = sorted(self.cations, key=lambda x: x.z) if self.there_are_cations() else None
        a = 0
        while (self.atoms.composition.fractional_composition != self.bulk.atoms.composition.fractional_composition and
               self.there_are_cations()):
            self.remove_sites([atomlist[a]])
            a = a * -1 if a < 0 else a * -1 - 1

def generate_mnx_slabs(bulk_str, hkl, thickness=15, vacuum=15, save=True):
    """
    Generates non-polar, stoichiometric slabs from a given
    bulk and set of miller indices. In addition, polyatomic
    anion clusters, which are bound covalently, are preserved.
    This is achieved by:

    1.  Generating a set of polar and non-polar surfaces using
        pymatgen's SlabGenerator class.
    2.  Attempts to Tasker-correct the slabs to achieve non-polarity,
        using pymatgen's get_tasker2_slabs function. This often finds
        more non-polar slabs.
    3.  Orthogonalises the slabs.
    4.  Corrects the slab to achieve non-polarity and bulk stoichiometry:
        i.      Removes lone anions which have been stripped from their atomic
                center.
        ii.     Since polarity is determined by the polyatomic anion distribution
                as it is the main non-divisible unit, a scaffold wherein all
                cations have been removed is generated.
        iii.    Non-polarity is achieved in the scaffold containing polyatomic
                anions.
        iv.     Non-polar scaffold is populated with the cations that
                were removed.
        v.      Non-polarity is achieved in the structure by selectively
                removing cations from one or the other side.
        vi.     Fractional composition of the slab is made to match that
                of the bulk by removing equivalent cations from both sides
                of the slab, while attempting to keep non-polarity.

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
    bulk_obj = BulkUnit(bulk)
    initial_slabs = tools.get_initial_slabs(bulk, hkl, thickness, vacuum)

    final_slabs = []
    for index, slab in enumerate(initial_slabs):
        load_bar(index, len(initial_slabs))
        initial_slab = slab.get_orthogonal_c_slab()
        scaffold = SlabUnit(initial_slab.copy(), bulk_obj)
        while len(np.append(scaffold.clusters_to_remove, scaffold.lone_anions)):
            scaffold.remove_sites(np.append(scaffold.clusters_to_remove, np.array(scaffold.lone_anions, dtype=object)))
            if center_str not in scaffold.atoms.formula:
                continue
        scaffold.remove_element(cations_strs)
        topbot = scaffold.depolarize_anions()
        reconstruction = SlabUnit(initial_slab.copy(), bulk_obj)
        reconstruction.remove_sites([site for site in reconstruction.atoms
                                        if site not in scaffold.atoms
                                        and site.element not in cations_strs])
        reconstruction.depolarize_cations(topbot)
        reconstruction.stoichiometrize()
        tools.set_site_attributes(reconstruction.atoms)
        if (not reconstruction.atoms.is_polar(tol_dipole_per_unit_area=0.01) and
            reconstruction.there_are_cations() and
            reconstruction.undercoordinated_sites(reconstruction.centers).size == 0):
            final_slabs.append(reconstruction)
        sys.stdout.flush()

    print('\nRemoving equivalent slabs...')
    tools.remove_equivalent_slabs(final_slabs)

    print(f'{len(final_slabs)} non-polar, stoichiometric slabs generated.')

    if save:
        for number, slab in enumerate(final_slabs):
            miller = f'{slab.atoms.miller_index[0]}{slab.atoms.miller_index[1]}{slab.atoms.miller_index[2]}'
            slab.atoms.to(f'{number}-{miller}.cif')

    return final_slabs

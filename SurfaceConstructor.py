from pymatgen.io.ase import AseAtomsAdaptor as ac
from pymatgen.core.surface import SlabGenerator, generate_all_slabs, get_symmetrically_distinct_miller_indices, center_slab
from pymatgen.analysis.structure_matcher import StructureMatcher
import pandas
import itertools
import os, copy, sys

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
        self.set_site_attributes(self.atoms)
        self.anion = self.anions[0].element
        self.center = self.centers[0].element
        self.cation = self.cations[0].element

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
        anions = [ atom for atom in self.atoms if atom.oxidation_state < 0 ]
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
        anion_neighbors = [ self.atoms.get_neighbors(atom, 2) for atom in self.anions ]
        anion_neighbors = [ neighbor[0] for neighbor in anion_neighbors if neighbor != [] ]
        center = list({object_.species_string: object_ for object_ in anion_neighbors}.values())
        if len(center) == 0:
            raise ValueError('No polyatomic anions could be found.')
        centers = [ atom for atom in self.atoms if atom.species_string == center[0].species_string ]
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
        all_cations = [ atom for atom in self.atoms
                            if atom.species_string != self.centers[0].species_string
                            and atom.oxidation_state > 0 ]
        if len(all_cations) == 0:
            raise ValueError('No cations could be found.')
        # unique_cations = list({object_.species_string: object_ for object_ in all_cations}.values())
        return all_cations

    @staticmethod
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
            distance_nn = min(
                                [structure.get_distance(site.index, atom.index)
                                 for atom in structure.get_neighbors(site, 3.5)]
                             ) + 0.2
            neighborlist = structure.get_neighbors(site, distance_nn)
            coordination_number = len(neighborlist)
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
            cluster.extend([ atom for atom in structure.get_neighbors(site, 2.3) ])
            cluster.append(site)
            return cluster

        #####
        ## Attributes of all pymatgen.core.sites.PeriodicSite objects
        ## conforming the structure are update here.
        ####
        for site in structure:
            site.index = structure.index(site)
            site.cluster = cluster(site, structure)
            site.coordination_number = coordination_number(site, structure)
            site.element = site.as_dict()['species'][0]['element']
            site.oxidation_state = site.as_dict()['species'][0]['oxidation_state']
        # for site in structure:
        #     max_coordination_number = max([atom.coordination_number for atom in structure if atom.element == site.element])
        #     site.undercoordinated = True if site.coordination_number < max_coordination_number else False

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
        self.set_site_attributes(self.atoms) # Update attributes of periodic sites.

    def coordination_numbers(self, sites):
        """
        Returns coordination numbers of a set of PeriodicSite objects
        in a list.

        Args:
            sites: list of PeriodicSite objects of which coordination numbers
            is to be obtained.

        Returns:
            (list): integers of coordination numbers.
        """
        return [ site.coordination_number for site in sites ]

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
        coord_no_bulk = [ site.coordination_number for site in self.bulk.atoms if site.element == sites[0].element ][0]
        # coord_no_bulk = len(self.bulk.atoms.get_neighbors(
        #                    [site for site in self.bulk.atoms if site.element == sites[0].element][0], 2.3)
        #                                                  )
        undercoordinated_sites = []
        for index, site in enumerate(sites):
            if site.coordination_number < coord_no_bulk:
                undercoordinated_sites.append(site)
        return undercoordinated_sites

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
        anions_list = [
                        atom for atom in self.anions
                             if self.bulk.centers[0].species_string not in
                             [site.species_string for site in atom.cluster]
                      ]
        return anions_list

    @property
    def clusters_to_remove(self):
        """
        Detects undercoordinated anionic centers and retrieves their nearest
        neighbors.

        Returns:
            (list): PeriodicSite objects of undercoordinated anionic centers,
            as well as their nearest anion neighbors.
        """
        atoms_to_remove = []
        for atom in self.undercoordinated_sites(self.centers):
            atoms_to_remove.extend(atom.cluster)
        return [ atom for atom in atoms_to_remove ]

    def remove_sites(self, sites):
        """
        Wrapping of the remove_sites function of pymatgen, this function
        removes a set of sites inside the SlabUnit class, without generating
        a copy of the structure and maintaining site attributes.

        Args:
            sites: list of PeriodicSite objects of sites which are to be
            removed from the main SlabUnit object.
        """
        template = self.atoms.copy()
        template.remove_sites([template.index(atom) for atom in template
                                                    if any(atom.is_periodic_image(site) for site in sites)])
        for site in [atom for atom in self.atoms if atom not in template]:
            self.atoms.remove(site)

    def remove_element(self, element):
        """
        Similarly to the 'remove_sites' function, this function removes
        all sites of a certain element.

        Args:
            element: element to be removed, as string (e.g., 'Fe').
        """
        atoms_list = [site for site in self.atoms if site.element == element]
        for atom in atoms_list:
            self.atoms.remove(atom)

    def top_site(self, element):
        """
        Finds the top site of a slab.

        Args:
            element: element of site to be investigated, as string (e.g., 'Fe').
        """
        atoms_list = [atom for atom in self.atoms if atom.element == element]
        atoms_list.sort(key=lambda x: x.z)
        return atoms_list[-1]

    def bottom_site(self, element):
        """
        Finds the bottom site of a slab.

        Args:
            element: element of site to be investigated, as string (e.g., 'Fe').
        """
        atoms_list = [atom for atom in self.atoms if atom.element == element]
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
               template.bulk.center in template.atoms.formula or
               len(template.centers)%2 != len(template.bulk.centers)%2 if template.bulk.center in template.atoms.formula
                                                                       else False):
                _sites_to_remove = template.top_site(template.bulk.center).cluster
                template.remove_sites(_sites_to_remove)
                sites_to_remove.extend(_sites_to_remove)
        if not template.bulk.center in template.atoms.formula:
            template = copy.deepcopy(self)
            sites_to_remove = []
            approach_value = 'bot'
            sites_to_remove = []
            while (template.atoms.is_polar(tol_dipole_per_unit_area=0.01) and
                   template.bulk.center in template.atoms.formula or
                   len(template.centers)%2 != len(template.bulk.centers)%2 if template.bulk.center in template.atoms.formula
                                                                           else False):
                _sites_to_remove = template.bottom_site(template.bulk.center).cluster
                template.remove_sites(_sites_to_remove)
                sites_to_remove.extend(_sites_to_remove)
        self.remove_sites(sites_to_remove) if template.atoms.formula != '' else None
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
        while (self.atoms.is_polar(tol_dipole_per_unit_area=0.01) and
               len(self.cations)/len(self.centers) > len(self.bulk.cations)/len(self.bulk.centers)):
            a += -1 if topbot == 'top' else 0
            self.remove_sites([atomlist[a]])
            a += 1 if topbot == 'bot' else 0

    def stoichiometrize(self):
        """
        Removes cations from a non-polar slab until composition
        matches the bulk stoichiometry.
        """
        atomlist = sorted(self.cations, key=lambda x: x.z) if self.bulk.cation in self.atoms.formula else None
        a = 0
        while (self.atoms.composition.fractional_composition != self.bulk.atoms.composition.fractional_composition and
               self.bulk.cation in self.atoms.formula):
            self.remove_sites([atomlist[a]])
            a = a * -1 if a < 0 else a * -1 - 1

def generate_slabs(bulk, hkl, thickness=15, vacuum=15):
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
    bulk.add_oxidation_state_by_guess()
    bulk_obj = BulkUnit(bulk)
    all_slabs = []
    print('Generating preliminary slabs...')
    if isinstance(hkl, int):
        hkl_list = get_symmetrically_distinct_miller_indices(bulk, hkl)
    if isinstance(hkl, list):
        hkl_list = hkl
    print(f'Miller indices that will be analysed: {", ".join([str(hkl) for hkl in hkl_list])}')
    for hkl in hkl_list:
        slabgen = SlabGenerator(bulk_obj.atoms, hkl, thickness, vacuum, center_slab=True)
        all_slabs.extend(slabgen.get_slabs())

    valid_slabs = []
    for slab in all_slabs:
        if not slab.is_polar(tol_dipole_per_unit_area=0.01):
            valid_slabs.append(slab)
        else:
            slab_tk_l = slab.get_tasker2_slabs(tol=0.01)
            for slab_tk in slab_tk_l:
                if not slab_tk.is_polar(tol_dipole_per_unit_area=0.01):
                    valid_slabs.append(slab_tk)
    for slab in all_slabs:
        if slab.is_polar(tol_dipole_per_unit_area=0.01):
            valid_slabs.append(slab)
    print(f'{len(valid_slabs)} preliminary slabs generated.')

    final_slabs = []
    for index, slab in enumerate(valid_slabs):
        percent = 100.*(index+1)/len(valid_slabs)
        sys.stdout.write('\r')
        sys.stdout.write("Sanitizing slabs: [{:{}}] {:>3}%"
                         .format('='*int(percent/(100.0/30)),
                                 30, int(percent)))
        initial_slab = slab.get_orthogonal_c_slab()

        scaffold = SlabUnit(initial_slab.copy(), bulk_obj)

        while (scaffold.clusters_to_remove + scaffold.lone_anions != [] if scaffold.bulk.center in scaffold.atoms.formula else False):
            scaffold.remove_sites(scaffold.clusters_to_remove + scaffold.lone_anions)
        if scaffold.bulk.center not in scaffold.atoms.formula:
            continue
        scaffold.remove_element(scaffold.bulk.cation)
        topbot = scaffold.depolarize_anions()
        reconstruction = SlabUnit(initial_slab.copy(), bulk_obj)
        reconstruction.remove_sites([site for site in reconstruction.atoms
                                          if site not in scaffold.atoms
                                          and site.element != reconstruction.bulk.cation])
        reconstruction.depolarize_cations(topbot)
        reconstruction.stoichiometrize()
        if (not reconstruction.atoms.is_polar(tol_dipole_per_unit_area=0.01) and
            reconstruction.bulk.cation in reconstruction.atoms.formula and
            reconstruction.undercoordinated_sites(reconstruction.centers) == []):
            reconstruction.set_site_attributes(reconstruction.atoms)
            final_slabs.append(reconstruction)
        sys.stdout.flush()

    print('\nRemoving equivalent slabs...')
    remove_equivalent_slabs(final_slabs)

    print(f'{len(final_slabs)} non-polar, stoichiometric slabs generated.')

    return final_slabs

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

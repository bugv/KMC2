"""Module containing the functions that manage structures and their representation

It contains functions to load structures, transform pymatgen structures to the corresponding arrays
and to modify these arrays.

Available functions:
- occupancy_vector_builder
- add_vacancy
- add_vacancy_random
- create_alloy
- atom_key_bulider
- empty_structure
- get_equiv_in_primative
- get_neighbour_from_displ
- alternate_neighbour_finder
- get_displ_one_site
- get_displacement_tensor
- get_max_nb_neighbours
- neighbour_finder (avoid this as order of neighbours is unknown)
- current_position_array_builder
- data_collector_builder
- time_colector_builder
- index_vector_builder
- find_vac
- get lattice_vectors
- get_comp_dict

Class: AtomPositions

Available class: AtomPositions
"""

import numpy as np
import pymatgen.core as pmg
from dataclasses import dataclass


def occupancy_vector_builder(structure: pmg.Structure, atom_key: dict) -> np.array:
    """Creates the occupancy vector where each index corresponds to a site in the structure and the value corresponds to the type of atom present

    :param structure: supercell
    :type structure: pmg.Structure
    :param atom_key: key linking type of atom to an int (use atom_key_builder)
    :type atom_key: dict
    :return: vector where values correspond to atom at the site with that index
    :rtype: np.array
    """
    occupancy_vector = np.empty(structure.num_sites)
    for i in range(structure.num_sites):
        occupancy_vector[i] = atom_key[
            str(structure[i].species.elements[0])
        ]  # index 0 because strcuture[i].species.elements is a list with only one entry
    return occupancy_vector


def add_vacancy(structure: pmg.Structure, vacancy_position: int) -> pmg.Structure:
    """Function to add a vacancy to the structure

    :param structure: supercell
    :type structure: pmg.Structure
    :param vacancy_position: position in structure to add vacancy
    :type vacancy_position: int
    :return: supercell with vacancy added
    :rtype: pmg.Structure
    """
    structure[vacancy_position] = pmg.DummySpecies("X")
    return structure


def add_vacancy_random(structure: pmg.Structure) -> tuple:
    """Function to randomly select the starting position of the vacancy

    :param structure: supercell
    :type structure: pmg.Structure
    :return: supercell with vacancy added and an int corresponding to the position of the vacancy
    :rtype: tuple
    """
    num_sites = structure.num_sites
    rng = np.random.default_rng()
    vacancy_position = int(np.floor(rng.random() * num_sites))
    structure[vacancy_position] = pmg.DummySpecies("X")
    return structure, vacancy_position


def create_alloy(structure: pmg.Structure, composition: dict) -> pmg.Structure:
    """Function to create a random alloy based on a composition

    :param structure: supercell in which to create the alloy
    :type structure: pmg.Structure
    :param composition: dict where the keys are the elements and the values their fractional composition
    :type composition: dict
    :raises ValueError: checks that the composition adds up to 1 (tolerance of 1e-4)
    :return: structure where types of atoms satisfy the requested composition
    :rtype: pmg.Structure
    """
    elem_array = np.array(list(composition.keys()))
    comp_array = np.array(list(composition.values()))
    if np.isclose(np.sum(comp_array), 1, rtol=1e-4, atol=1e-3) == False:
        raise ValueError(
            "The sum of the composition of the alloy is not 1. Check the composition of the alloy"
        )
    num_sites_elem = np.floor(comp_array * structure.num_sites).astype(int)
    if np.sum(num_sites_elem) != structure.num_sites:
        num_atoms_to_add = structure.num_sites - np.sum(num_sites_elem)
        np.argmax(num_sites_elem)
        num_sites_elem[np.argmax(num_sites_elem)] += num_atoms_to_add
    list_sites = np.arange(structure.num_sites)
    random_order_sites = np.random.permutation(list_sites)
    elements_repeated = np.repeat(elem_array, num_sites_elem)
    for site in range(structure.num_sites):
        structure[int(random_order_sites[site])] = elements_repeated[site]
    for elem in elem_array:
        print("elem test", elem)
        print(
            "Composition of the structure created:",
            elem,
            " = ",
            structure.composition.get_atomic_fraction(elem),
        )
    return structure


def atom_key_builder(structure: pmg.Structure) -> dict:
    """Create a dict to link atom type to an integer

    :param structure: supercell
    :type structure: pmg.Structure
    :return: dict, with types of atoms as key and unique ints as values
    :rtype: dict
    """
    list_values = list(range(len(structure.species)))
    element_list = structure.elements
    for i in range(
        len(element_list)
    ):  # NOTE: not ideal to use range of len but I couldn't enumerate over elements in list to convert them to string
        element_list[i] = str(element_list[i])
    atom_key = dict(zip(element_list, list_values))
    return atom_key


def empty_structure(supercell: pmg.Structure) -> pmg.Structure:
    """Creates a supercell where all sites are occupied by vacancies

    :param supercell: supercell
    :type supercell: pmg.Structure
    :return: same supercell where all atoms are replaced by vacancies
    :rtype: pmg.Structure
    """
    for site in supercell:
        site.species = pmg.DummySpecies("X")
    return supercell


def get_equiv_in_primative(
    primitivecell: pmg.Structure, supercell: pmg.Structure
) -> np.array:
    """Create an array which indicates which site in the primitive cell,
       the site in the supercell is equivalent to

    :param supercell: supercell
    :type supercell: pmg.Structure
    :raises ValueError: check that a site is not equivlent to two sites in the primitive cell
    :raises ValueError: check that all sites in the supercell are equivalent to a site in the primitive cell
    :return: array where at each index there is the equivalent site in the primitive cell
    :rtype: np.array
    """

    lattice_matrix_inv = np.linalg.inv(np.transpose(primitivecell.lattice.matrix))
    equivalent_sites_array = np.full(supercell.num_sites, None)
    for i in range(supercell.num_sites):
        for j in range(primitivecell.num_sites):
            result = lattice_matrix_inv @ (
                supercell.cart_coords[i, :] - primitivecell.cart_coords[j, :]
            )
            is_equiv = np.all(np.isclose(result, np.round(result), atol=1e-5))
            if is_equiv and equivalent_sites_array[i] is not None:
                raise ValueError(
                    "site is equivalent to two sites in the primitive cell"
                )
            if equivalent_sites_array[i] is None and is_equiv:
                equivalent_sites_array[i] = j
    return equivalent_sites_array


def get_neighbours_from_displ(
    structure: pmg.Structure,
    equivalent_sites_array: np.array,
    radius: float,
    disp_tensor: np.array,
) -> np.array:
    """Creates a neighbour array from the displacement vectors and the primitive equivalent sites

    :param structure: supercell
    :type structure: pmg.Structure
    :param equivalent_sites_array: array with equivelent sites in primitive cell
    :type equivalent_sites_array: np.array
    :param radius: cutoff radius for neighbours
    :type radius: float
    :param disp_tensor: rank three tensor with displacements
    :type disp_tensor: np.array
    :raises ValueError: neighbour site from displacement not found in the supercell
    :return: neighbour array, where each column contains the neighbours of the atom at that index
    :rtype: np.array
    """
    structure = empty_structure(structure)
    neighbours_array = np.full(
        (get_max_nb_neighbours(structure, radius), structure.num_sites), None
    )
    coords_array = structure.cart_coords
    sites_list = structure.sites
    for center_index, center_site in enumerate(coords_array):
        prim_equiv = equivalent_sites_array[center_index]
        displacements = disp_tensor[:, :, prim_equiv]
        for nb_neighbour, displ in enumerate(np.transpose(displacements)):
            if np.all(displ != None):
                new_frac_coords = np.mod(
                    np.round(
                        np.array(
                            structure.lattice.get_fractional_coords(
                                center_site + displ
                            ),
                            dtype=float,
                        ),
                        6,
                    ),
                    1,
                )
                new_site = pmg.PeriodicSite(
                    species="X",
                    coords=np.mod(new_frac_coords, 1),
                    lattice=structure.lattice,
                    to_unit_cell=True,
                    coords_are_cartesian=False,
                )
                index = -1
                for site in structure:
                    test_dist = site.distance(new_site) <= 1e-4
                    if test_dist:
                        index = sites_list.index(site)
                        break
                if index == -1:
                    raise ValueError("neighbour not found with alternative method")
                neighbours_array[nb_neighbour, center_index] = index
    return neighbours_array


def alternate_neighbour_finder(supercell: pmg.Structure, radius: float) -> np.array:
    """Function to create an array containing the neighbours of a site. Method based on always ordering the neighbours in the same order based on their displacements relative to the central site.

    :param supercell: supercell
    :type supercell: pmg.Structure
    :param radius: cutoff radius for neighbours
    :type radius: float
    :return: Array where each columns contains the index of the neighbours of the site with the same index as the column.
    :rtype: np.array
    """
    max_nb_neighbours = get_max_nb_neighbours(supercell, radius)
    nb_sites = supercell.num_sites
    neighbour_array = np.full((max_nb_neighbours, nb_sites), None)
    csites, nsites, offsets, dists = supercell.get_neighbor_list(radius)
    cart_nsites = np.array([supercell.sites[i].coords for i in nsites])
    cart_csites = np.array([supercell.sites[i].coords for i in csites])
    lattice = get_lattice_vectors(supercell).transpose()
    displs = cart_nsites - cart_csites + offsets @ lattice
    sorter = np.lexsort(
        (
            displs[:, 2],
            displs[:, 1],
            displs[:, 0],
        )
    )
    sorted_displs = displs[sorter]
    sorted_csites = csites[sorter]
    sorted_nsites = nsites[sorter]
    for index in range(np.shape(neighbour_array)[1]):
        neighbour_array[:, index] = sorted_nsites[np.where(sorted_csites == index)]
    return neighbour_array


def get_displ_one_site(site: int, structure: pmg.Structure, radius: float) -> tuple:
    """Function that gets the displacement vectors for a given site in a structure

    :param site: index of the site in the structure
    :type site: int
    :param structure: supercell
    :type structure: pmg.Structure
    :param radius: cutoff radius for neighbours
    :type radius: float
    :return: 3x max nb neighbours in the structure array where each column is a displacement vector
             + array with indices of neighbours in the same order
    :rtype: tuple
    """
    displacements = np.full((3, get_max_nb_neighbours(structure, radius)), None)
    neighbour_list = structure.get_neighbor_list(radius)
    lattice = structure.lattice.matrix
    coords_site = structure[site].coords
    center_indices = np.where(neighbour_list[0] == site)
    indices_neighbours = np.full((1, get_max_nb_neighbours(structure, radius)), None)
    for index_to_add, index_in_list in np.ndenumerate(center_indices):
        offset = neighbour_list[2][index_in_list]
        neighbour_coords = (
            structure[neighbour_list[1][index_in_list]].coords + offset @ lattice
        )
        displ = neighbour_coords - coords_site
        displacements[:, index_to_add[1]] = displ
        indices_neighbours[0, index_to_add[1]] = neighbour_list[1][index_to_add[1]]
    return (displacements, indices_neighbours)


def get_displacement_tensor(
    primitive: pmg.Structure, supercell: pmg.Structure, radius: float
) -> np.array:
    """Create a rank three tensor of size 3 x max nb of neighbour x nb sites primitive cell
       with the displacement for each site in the primitive cell.
       Where each layer contains the displacement vectors for one site in the primitive cell,
       and each column is a diplacement vector

    :param primitive: primitive cell
    :type primitive: pmg.Structure
    :param supercell: supercell
    :type supercell: pmg.Structure
    :param radius: cutoff radius ffor neighbours
    :type radius: float
    :return: randk three tensor of size   max nb neighbours x 3 x nb sites primitive cell,
    where each layer contains the displacement vectors for one site in the primitive cell,
    and each row is a displacement vector.
    :rtype: np.array
    """
    neighbour_list = primitive.get_neighbor_list(radius)
    displ_vect_tensor = np.full(
        (
            3,
            get_max_nb_neighbours(supercell, radius),
            primitive.num_sites,
        ),
        None,
    )
    index_list = np.array([])
    for i in range(primitive.num_sites):

        count = np.count_nonzero(neighbour_list[0] == i)
        indices = np.arange(0, count, 1)
        index_list = np.append(index_list, indices)
    lattice = primitive.lattice.matrix
    for count_center, center in enumerate(neighbour_list[0]):
        offset = neighbour_list[2][count_center]
        neighbour_coords = primitive[neighbour_list[1][count_center]].coords
        neighbour_coords = neighbour_coords + offset @ lattice
        displ = neighbour_coords - primitive[center].coords
        displ_vect_tensor[:, int(index_list[count_center]), center] = displ
    displ_vect_tensor = np.transpose(displ_vect_tensor, (1, 0, 2))
    sorted_displ_tensor = np.empty_like(displ_vect_tensor)
    for i in range(np.shape(sorted_displ_tensor)[2]):
        layer = displ_vect_tensor[:, :, i]
        sorter = np.lexsort((layer[:, 2], layer[:, 1], layer[:, 0]))
        sorted_displ_tensor[:, :, i] = layer[sorter]
    sorted_displ_tensor = np.transpose(sorted_displ_tensor, (1, 0, 2))
    return sorted_displ_tensor


def get_max_nb_neighbours(structure: pmg.Structure, radius: int) -> int:
    """get the maximum number of neighbours in the structure

    :param structure: supercell
    :type structure: pmg.Structure
    :param radius: cutoff radius for neighbours
    :type radius:
    :return: the maximum number of neighbours in the structure
    :rtype: int
    """
    neighbour_list = structure.get_neighbor_list(radius)
    centers_list = neighbour_list[0]
    nb_sites_struct = structure.num_sites
    counts = np.full(nb_sites_struct, None)
    for i in range(nb_sites_struct):
        counts[i] = np.count_nonzero(centers_list == i)
    max_neighbours = np.max(counts)
    return max_neighbours


def neighbour_finder(structure: pmg.Structure, radius: float) -> np.array:
    """A function to find the neighbouring sites for each site

    :param structure: A supercell
    :type structure: pmg.Structure
    :param radius: THe cutoff radius for a neighbour
    :type radius: float
    :return: An array where the values in a given column contains the labels for the neighbouring sites
    :rtype: np.array
    """
    N = structure.num_sites
    neighbour_list = structure.get_neighbor_list(radius)
    centers_list = neighbour_list[0]
    nearest_sites_list = neighbour_list[1]

    M = get_max_nb_neighbours(structure, radius)
    neighbour_array = np.full((M, N), None)
    for index, atom in enumerate(centers_list):
        added_to_list = False
        for m in range(M):
            if neighbour_array[m, atom] == None and not added_to_list:
                neighbour_array[m, atom] = nearest_sites_list[index]
                added_to_list = True
    return neighbour_array


def current_position_array_builder(supercell: pmg.Structure) -> np.array:
    """Create an array to store the position of the atoms at the current time

    :param supercell: supercell
    :type supercell: pmg.Structure
    :return: An array with the coordinates of the atoms at the initial time
    :rtype: np.array
    """
    return supercell.cart_coords.transpose()


def data_collector_builder(
    supercell: pmg.Structure, sampling_freq: float, total_number_step: int
) -> np.array:
    """Create an array in which to store the position of the atoms at the sampling times
    Each layer is a sampling step (in time), and in each layer, each column corresponds
    to the coordinates of the atom that started in that site
    :param supercell: supercell
    :type supercell: pmg.Structure
    :param sampling_freq: Frequency (nb steps) for sampling atom position
    :type sampling_freq: float
    :param total_number_step: Number of steps
    :type total_number_step: int
    :return: 3D array to store the results at each sampling step
    :rtype: np.array
    """
    data_collector = np.full(
        (3, supercell.num_sites, int(total_number_step / sampling_freq)), None
    )
    data_collector[:, :, 0] = supercell.cart_coords.transpose()
    return data_collector


def time_collector_builder(sampling_freq: float, total_number_step: int) -> np.array:
    """Function to create an array to store the times for sampling

    :param sampling_freq: frequency at which to take samples (nb of steps)
    :type sampling_freq: float
    :param total_number_step: Total number of iterations of the KMC
    :type total_number_step: int
    :return: array to store the times at which the structure is sampled
    :rtype: np.array
    """
    return np.full((1, int(total_number_step / sampling_freq)), 0.0)


def index_vector_builder(structure: pmg.Structure) -> np.array:
    """Initialize an index vector which links the occupancy vector(species at each site)
       and the current position of the atoms based on their starting site

    :param structure: supercell
    :type structure: pmg.Structure
    :return: 1D array with position of atoms at the initial situation
    :rtype: np.array
    """
    return np.arange(structure.num_sites)


def find_vac(occupancy_vector: np.array, atom_key: dict) -> int:
    """Funtion to find the position of the vacancy in the occupancy vector

    :param occupancy_vector: occupancy vector
    :type occupancy_vector: np.array
    :param atom_key: atom key
    :type atom_key: dict
    :return: integer correponding to the index of the vacancy within the occupancy vector
    :rtype: int
    """
    vac_pos = np.where(occupancy_vector == atom_key["X0+"])
    return int(vac_pos[0])


def get_lattice_vectors(structure: pmg.Structure) -> np.array:
    """Function to get the lattice vectors forprint(csites)
    print(cart_csites) the supercell

    :param structure: supercell
    :type structure: pmg.Structure
    :return: 3 x 3 array with each column a lattice vector
    :rtype: np.array
    """
    lattice = structure.lattice.as_dict()
    return np.vstack(lattice["matrix"]).transpose()


def get_comp_dict(supercell: pmg.Structure, atomkey: dict) -> dict:
    """Function to create a dict containing the composition of the supercell

    :param supercell: supercell
    :type supercell: pmg.Structure
    :param atomkey: atom key dict created by the atom_key_builder function, should correspond to the same supercell
    :type atomkey: dict
    :return: Dict where the keys are the elements and the values are the atomic fraction of the element in the supercell
    :rtype: dict
    """
    element_list = list(atomkey.keys())
    comp_dict = {}
    for element in element_list:
        comp_dict[element] = supercell.composition.get_atomic_fraction(element)
    return comp_dict


@dataclass
class AtomPositions:
    """Class to manage all the arrays that describe the position of the atoms
    Consists of occupancy_vector, current_position_array, index_array
    """

    occupancy_vector: np.array
    current_position_array: np.array
    index_array: np.array
    neighbour_list: np.array
    equivalent_sites_in_prim: np.array
    displacement_tensor: np.array

    def swap(self, i: int, j: int):
        """Method to swap the species in sites a and b of the structure

        :param index_a: index of the first species
        :type index_a: int
        :param index_b: index of the second species
        :type index_b: int
        """

        self.occupancy_vector[[i, j]] = self.occupancy_vector[[j, i]]
        a = np.where(self.index_array == i)
        b = np.where(self.index_array == j)
        temp_a = np.copy(self.current_position_array[:, a])
        self.current_position_array[:, a] = self.current_position_array[:, b]
        self.current_position_array[:, b] = temp_a
        temp_i = self.index_array[i]
        self.index_array[i] = self.index_array[j]
        self.index_array[j] = temp_i

    def swap_displ(
        self,
        initial: int,
        final: int,
    ):
        neighbours_of_initial = self.neighbour_list[:, initial]
        position_final_in_list = np.where(neighbours_of_initial == final)
        if np.size(position_final_in_list) == 0:
            raise ValueError("The two sites are not neighbours")
        prim_equiv = self.equivalent_sites_in_prim[initial]
        initial_start_index = self.index_array[initial]
        final_start_index = self.index_array[final]
        displ = self.displacement_tensor[:, position_final_in_list, prim_equiv]
        self.current_position_array[0, initial_start_index] = (
            self.current_position_array[0, initial_start_index] + displ[0]
        )
        self.current_position_array[1, initial_start_index] = (
            self.current_position_array[1, initial_start_index] + displ[1]
        )
        self.current_position_array[2, initial_start_index] = (
            self.current_position_array[2, initial_start_index] + displ[2]
        )
        self.current_position_array[0, final_start_index] = (
            self.current_position_array[0, final_start_index] - displ[0]
        )
        self.current_position_array[1, final_start_index] = (
            self.current_position_array[1, final_start_index] - displ[1]
        )
        self.current_position_array[2, final_start_index] = (
            self.current_position_array[2, final_start_index] - displ[2]
        )
        self.occupancy_vector[[initial, final]] = self.occupancy_vector[
            [final, initial]
        ]

        self.index_array[[initial, final]] = self.index_array[[final, initial]]

    def swap_displ_V2(
        self,
        initial: int,
        position_final_in_list: int,
    ):
        neighbours_of_initial = self.neighbour_list[:, initial]
        final = neighbours_of_initial[position_final_in_list]
        prim_equiv = self.equivalent_sites_in_prim[initial]
        initial_start_index = self.index_array[initial]
        final_start_index = self.index_array[final]
        displ = self.displacement_tensor[:, position_final_in_list, prim_equiv]
        self.current_position_array[0, initial_start_index] = (
            self.current_position_array[0, initial_start_index] + displ[0]
        )
        self.current_position_array[1, initial_start_index] = (
            self.current_position_array[1, initial_start_index] + displ[1]
        )
        self.current_position_array[2, initial_start_index] = (
            self.current_position_array[2, initial_start_index] + displ[2]
        )
        self.current_position_array[0, final_start_index] = (
            self.current_position_array[0, final_start_index] - displ[0]
        )
        self.current_position_array[1, final_start_index] = (
            self.current_position_array[1, final_start_index] - displ[1]
        )
        self.current_position_array[2, final_start_index] = (
            self.current_position_array[2, final_start_index] - displ[2]
        )
        self.occupancy_vector[[initial, final]] = self.occupancy_vector[
            [final, initial]
        ]
        self.index_array[[initial, final]] = self.index_array[[final, initial]]

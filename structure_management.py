"""Module containing the functions that manage structures and their representation

It contains functions to load structures, transform pymatgen structures to the corresponding arrays
and to modify these arrays.

Available functions:
- occupancy_vector_builder
- atom_key_bulider
- empty_structure
- get_equiv_in_primative
- get_neighbour_from_displ
- get_displacement_tensor
- get_max_nb_neighbours
- neighbour_finder (avoid this as order of neighbours is unknown)
- current_position_array_builder
- data_collector_builder
- time_colector_builder
- index_vector_builder
- find_vac

Available class: AtomPositions
"""

import numpy as np
import pymatgen.core as pmg
from dataclasses import dataclass
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import psutil


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
    for i in range(structure.num_sites):  # can't enumerate over empty array
        occupancy_vector[i] = atom_key[
            str(structure[i].species.elements[0])
        ]  # index 0 because strcuture[i].species.elements is a list with only one entry
    return occupancy_vector


def atom_key_builder(structure: pmg.Structure) -> dict:
    """Create a dict to link atom type to an integer

    :param structure: supercell
    :type structure: pmg.Structure
    :return: dict, with types of atoms as key and unique ints as values
    :rtype: dict
    """
    list_values = list(range(len(structure.species)))
    element_list = structure.species
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


def get_equiv_in_primative(supercell: pmg.Structure) -> np.array:
    """Create an array which indicates which site in the primitive cell,
       the site in the supercell is equivalent to

    :param supercell: supercell
    :type supercell: pmg.Structure
    :raises ValueError: check that a site is not equivlent to two sites in the primitive cell
    :raises ValueError: check that all sites in the supercell are equivalent to a site in the primitive cell
    :return: array where at each index there is the equivalent site in the primitive cell
    :rtype: np.array
    """
    supercell = empty_structure(supercell)
    primitivecell = SpacegroupAnalyzer(supercell).find_primitive()
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
    # if np.any(equivalent_sites_array, None) != 1:
    #     print("check equivalent sites array", equivalent_sites_array)
    #     print("test condition", np.any(equivalent_sites_array, None))
    #     raise ValueError("one or more sites has no equivalent in the primitive cell")
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
        process = psutil.Process()
        mem_before = process.memory_info().rss / 1024  # in kilobytes
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
                if index == -1:
                    raise ValueError("neighbour not found with alternative method")
                neighbours_array[nb_neighbour, center_index] = index
        mem_after = process.memory_info().rss / 1024  # in kilobytes
        mem_used = mem_after - mem_before
        print("mem used for one loop neighbour finder, all neighbours one atom (kilobytes)", mem_used)

    return neighbours_array


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


def get_displacement_tensor(structure: pmg.Structure, radius: float) -> np.array:
    """Create a rank three tensor of size 3 x max nb of neighbour x nb sites primitive cell
       with the displacement for each site in the primitive cell.
       Where each layer contains the displacement vectors for one site in the primitive cell,
       and each column is a diplacement vector

    :param structure: supercell
    :type structure: pmg.Structure
    :param radius: cutoff radius ffor neighbours
    :type radius: float
    :return: randk three tensor of size 3 x max nb neighbours x nb sites primitive cell,
    where each layer contains the displacement vectors for one site in the primitive cell,
    and each column is a displacement vector.
    :rtype: np.array
    """
    structure = empty_structure(structure)
    primitive = SpacegroupAnalyzer(structure).find_primitive()
    neighbour_list = primitive.get_neighbor_list(radius)
    displ_vect_tensor = np.full(
        (
            3,
            get_max_nb_neighbours(structure, radius),
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
    return displ_vect_tensor


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


def current_position_array_builder(structure: pmg.Structure) -> np.array:
    """Create an array to store the position of the atoms at the current time

    :param structure: supercell
    :type structure: pmg.Structure
    :return: An array with the coordinates of the atoms at the initial time
    :rtype: np.array
    """
    return structure.cart_coords.transpose()


def data_collector_builder(
    structure: pmg.Structure, sampling_freq: float, total_number_step: int
) -> np.array:
    """Create an array in which to store the position of the atoms at the sampling times
    Each layer is a sampling step (in time), and in each layer, each column corresponds
    to the coordinates of the atom that started in that site
    :param structure: supercell
    :type structure: pmg.Structure
    :param sampling_freq: Frequency (nb steps) for sampling atom position
    :type sampling_freq: float
    :param total_number_step: Number of steps
    :type total_number_step: int
    :return: 3D array to store the results at each sampling step
    :rtype: np.array
    """
    data_collector = np.full(
        (3, structure.num_sites, int(total_number_step / sampling_freq)), None
    )
    data_collector[:, :, 0] = structure.cart_coords.transpose()
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
    """Function to get the lattice vectors for the supercell

    :param structure: supercell
    :type structure: pmg.Structure
    :return: 3 x 3 array with each column a lattice vector
    :rtype: np.array
    """
    lattice = structure.lattice.as_dict()
    return np.vstack(lattice["matrix"]).transpose()


@dataclass
class AtomPositions:
    """Class to manage all the arrays that describe the position of the atoms
    Consists of occupancy_vector, current_position_array, index_array
    """

    occupancy_vector: np.array
    current_position_array: np.array
    index_array: np.array
    # Unsure to include the following
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

        # input_seq[[ix1, ix2]] = input_seq[[ix2, ix1]]

        # update occupancy vector
        self.occupancy_vector[[i, j]] = self.occupancy_vector[[j, i]]
        # temp_i = self.occupancy_vector[i]
        # self.occupancy_vector[i] = self.occupancy_vector[j]
        # self.occupancy_vector[j] = temp_i

        # update current_position_array
        a = np.where(self.index_array == i)
        b = np.where(self.index_array == j)
        temp_a = np.copy(self.current_position_array[:, a])
        self.current_position_array[:, a] = self.current_position_array[:, b]
        self.current_position_array[:, b] = temp_a

        # update index array
        temp_i = self.index_array[i]
        self.index_array[i] = self.index_array[j]
        self.index_array[j] = temp_i

    # New  swap including the displacement vectors
    # Input: indices of the two atoms/things to swap, list of equivelanet sites in primitive,
    #        displacement vector tensor
    # Output: nothing, modifies the atom position object
    # update the positions using the displacements vectors
    # for the first site, look through the neighbour list to get which vertical index the second site is at
    # Get wich primitive site the first site is equivalent to
    # get the displacement vector
    # add the diplacemnt vector to the current position of the atom1
    # subtract the displacement vector to the current position of atom 2
    # update the occupancy vector (swap the )
    # update the index array
    ## Question: should the atoms position object contain the list of equivalent sites in the primitive
    # and the displacement vectors and the neighbours list
    # Pro: call swap only acts on self (no need to input displacements...),
    #       lower risk of copying more than necessary
    # Con: Not really part of the atoms positions

    # Should it contain a check that the two sites are actually neighbours
    # Has to contain it because update position is based on them being neighbours
    def swap_displ(
        self,
        initial: int,
        final: int,
    ):
        # print("sites swaped", initial, final)
        # get the neighbour list for the initial site
        neighbours_of_initial = self.neighbour_list[:, initial]
        # find where in the neighbour list the final site is
        # np.where not ideal, but looking for an int in a small number
        position_final_in_list = np.where(neighbours_of_initial == final)
        # Check site is actually a neighbour
        if np.size(position_final_in_list) == 0:
            raise ValueError("The two sites are not neighbours")
        # Get the prim equiv for initial site
        prim_equiv = self.equivalent_sites_in_prim[initial]
        # get where the initial is in the index array
        # TODO is this possible without using the np.where which is a for loop behind the scenes
        initial_start_index = self.index_array[initial]
        # print("pos_inital (should always be the same)", initial_start_index)
        # get where the final is in the index array
        final_start_index = self.index_array[final]
        # use position in neighour list and prim equiv to get displ between initial and final
        displ = self.displacement_tensor[:, position_final_in_list, prim_equiv]
        # print("displ", displ)
        # add displ to position of initial in current position array
        # for first coord
        self.current_position_array[0, initial_start_index] = (
            self.current_position_array[0, initial_start_index] + displ[0]
        )
        # for second coord
        self.current_position_array[1, initial_start_index] = (
            self.current_position_array[1, initial_start_index] + displ[1]
        )
        # for third coord
        self.current_position_array[2, initial_start_index] = (
            self.current_position_array[2, initial_start_index] + displ[2]
        )
        # update position of final by subtracting displ from coords
        # for first coord
        self.current_position_array[0, final_start_index] = (
            self.current_position_array[0, final_start_index] - displ[0]
        )
        # for second coord
        self.current_position_array[1, final_start_index] = (
            self.current_position_array[1, final_start_index] - displ[1]
        )
        # for third coord
        self.current_position_array[2, final_start_index] = (
            self.current_position_array[2, final_start_index] - displ[2]
        )
        # swap in occupancy vector
        self.occupancy_vector[[initial, final]] = self.occupancy_vector[
            [final, initial]
        ]

        # update index array by swapping the index array values
        self.index_array[[initial, final]] = self.index_array[[final, initial]]

        # neighbour_position_final_rel_to_initial = n  # What is this line doing

    # def swap2(self, i: int, j: int):
    #     """Method to swap the species in sites a and b of the structure

    #     :param index_a: index of the first species
    #     :type index_a: int
    #     :param index_b: index of the second species
    #     :type index_b: int
    #     """
    #     # update occupancy vector
    #     temp_i = self.occupancy_vector[i]
    #     self.occupancy_vector[i] = self.occupancy_vector[j]
    #     self.occupancy_vector[j] = temp_i

    #     # update current_position_array
    #     a = np.where(self.index_array == i)
    #     b = np.where(self.index_array == j)

    #     distance = self.frac_coord_array[:a] - self.frac_coord_array[:b]
    #     a_displacement = [0, 0, 0]
    #     b_displacement = [0, 0, 0]
    #     for k in range(3):
    #         if abs(distance[k]) > 0.5:
    #             a_displacement = (1 - abs(distance[k])) * self.lattice_vectors[:k]
    #             b_displacement = -(1 - abs(distance[k])) * self.lattice_vectors[:k]

    #     self.current_position_array[:, a] = (
    #         self.current_position_array[:, a] + a_displacement
    #     )
    #     self.current_position_array[:, b] = (
    #         self.current_position_array[:, b] + b_displacement
    #     )

    #     # update index array
    #     temp_i = self.index_array[i]
    #     self.index_array[i] = self.index_array[j]
    #     self.index_array[j] = temp_i

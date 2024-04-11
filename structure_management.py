"""Module containing the functions that manage structures and their representation

It contains functions to load structures, transform pymatgen structures to the corresponding arrays
and to modify these arrays.

Available functions:
- atom_key_bulider
- occupancy_vector_builder
- neighbour_finder"""

import numpy as np
import pymatgen.core as pmg
from dataclasses import dataclass
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

## Read from file to structure
# FUnction to read structure from file and produce a pymatgen structure
# Input: a file location
# Output: a pymatgen structure


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
    list_values = list(range(len(structure.elements)))
    element_list = structure.elements
    for i in range(
        len(element_list)
    ):  # NOTE: not ideal to use range of len but I couldn't enumerate over elements in list to convert them to string
        element_list[i] = str(element_list[i])
    atom_key = dict(zip(element_list, list_values))
    return atom_key


def get_equiv_primative(supercell: pmg.Structure) -> np.array:
    """Create an array which indicates which site in the primitive cell,
       the site in the supercell is equivalent to

    :param supercell: supercell
    :type supercell: pmg.Structure
    :raises TypeError: check that a site is not equivlent to two sites in the primitive cell
    :raises TypeError: check that all sites in the supercell are equivalent to a site in the primitive cell
    :return: array where at each index there is the equivalent site in the primitive cell
    :rtype: np.array
    """
    for site in supercell:
        site.species = pmg.DummySpecies("X")
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
                raise TypeError("site is equivalent to two sites in the primitive cell")
            if equivalent_sites_array[i] is None and is_equiv:
                equivalent_sites_array[i] = j
    if np.any(equivalent_sites_array, None) is 0.0:
        raise TypeError("one or more sites has no equivalent in the primitive cell")
    return equivalent_sites_array


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

    M = 0
    for index, atom in enumerate(centers_list):
        if centers_list[index] == 0:
            M += 1
    neighbour_array = np.full((M, N), None)
    for index, atom in enumerate(centers_list):
        added_to_list = False
        for m in range(M):
            if neighbour_array[m, atom] == None and not added_to_list:
                neighbour_array[m, atom] = nearest_sites_list[index]
                added_to_list = True
    return neighbour_array


def neighbour_finder_v2(structure: pmg.Structure, radius: float) -> np.array:
    neighbour_list = structure.get_neighbor_list(radius)
    centers_list = neighbour_list[0]
    nearest_sites_list = neighbour_list[1]
    unique, counts = np.unique(centers_list, return_counts=True)
    m = counts.max()
    neighbour_array = np.full((m, structure.num_sites), None)  # create neighbour array
    position_array = np.tile(np.arange(m), (structure.num_sites, 1)).T
    position_array[:, :] = np.where(
        np.arange(m) < counts[:, None], np.arange(m), None
    ).T
    # position_array = np.full((m, structure.num_sites), None)
    # position_array[:, :] = np.arange(0, m)
    print(position_array)
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


def frac_coord_array_builder(structure: pmg.Structure) -> np.array:
    """Function to get fractional coordinates of each site in structure

    :param structure: supercell
    :type structure: pmg.Structure
    :return: 3 x nb_sites array with fractional coordinates for each site,
    :rtype: np.array
    """
    return structure.frac_coords.transpose()


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
    frac_coord_array: dict
    lattice_vectors: np.array

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

    def swap2(self, i: int, j: int):
        """Method to swap the species in sites a and b of the structure

        :param index_a: index of the first species
        :type index_a: int
        :param index_b: index of the second species
        :type index_b: int
        """
        # update occupancy vector
        temp_i = self.occupancy_vector[i]
        self.occupancy_vector[i] = self.occupancy_vector[j]
        self.occupancy_vector[j] = temp_i

        # update current_position_array
        a = np.where(self.index_array == i)
        b = np.where(self.index_array == j)

        distance = self.frac_coord_array[:a] - self.frac_coord_array[:b]
        a_displacement = [0, 0, 0]
        b_displacement = [0, 0, 0]
        for k in range(3):
            if abs(distance[k]) > 0.5:
                a_displacement = (1 - abs(distance[k])) * self.lattice_vectors[:k]
                b_displacement = -(1 - abs(distance[k])) * self.lattice_vectors[:k]

        self.current_position_array[:, a] = (
            self.current_position_array[:, a] + a_displacement
        )
        self.current_position_array[:, b] = (
            self.current_position_array[:, b] + b_displacement
        )

        # update index array
        temp_i = self.index_array[i]
        self.index_array[i] = self.index_array[j]
        self.index_array[j] = temp_i

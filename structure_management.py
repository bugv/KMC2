"""Module containing the functions that manage structures and their representation

It contains functions to load structures, transform pymatgen structures to the corresponding arrays
and to modify these arrays.

Available functions:
- atom_key_bulider
- occupancy_vector_builder
- neighbour_finder"""

import numpy as np
import pymatgen.core as pmg

## Read from file to structure
# FUnction to read structure from file and produce a pymatgen structure
# Input: a file location
# Output: a pymatgen structure

## Create the occupancy vector
# Input: structure
# Output: Occupation vector of the structure with one vacancy in the structure (size=nb of sites in the structure)
# 1. read structure from the supercell
# 4. Get number of sites in the structure
# 5. Create a vector with length number of sites in the structure
# 6. Convert structure to occupancy vector by reading each item in the structure, and using key to get appropriate value to write to the vector
# 7. Return occupancy vector


def occupancy_vector_builder(structure: pmg.Structure, atom_key: dict) -> np.array:
    """Creates the occupancy vector where each index corresponds to a site in the structure and the value corresponds to the type of atom present

    :param structure: supercell
    :type structure: pmg.Structure
    :param atom_key: key linking type of atom to an int (use atom_key_builder)
    :type atom_key: dict
    :return: vector where values correspond to atom at the site with that index
    :rtype: np.array
    """
    occupancy_vector = np.empty(
        structure.num_sites
    )  # create empty occupancy vector of length  = nb of sites
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


def neighbour_finder(structure: pmg.Structure, radius: float) -> np.array:
    """A function to find the neighbouring sites for each site

    :param structure: A supercell
    :type structure: pmg.Structure
    :param radius: THe cutoff radius for a neighbour
    :type radius: float
    :return: An array where the values in a given column contains the labels for the neighbouring sites
    :rtype: np.array
    """
    N = structure.num_sites  # Get number of sites N
    neighbour_list = structure.get_neighbor_list(
        radius
    )  # NOTE cutoff value based on radius (but the value of radius is arbitrary, I used 2.5)
    centers_list = neighbour_list[0]
    nearest_sites_list = neighbour_list[1]
    # Get number of neighbours per site based on the number of neighbours of the first site
    M = 0
    for index, atom in enumerate(centers_list):
        if centers_list[index] == 0:
            M += 1
    # Create neighbours array and fill it
    neighbour_array = np.full((M, N), None)
    for index, atom in enumerate(centers_list):
        added_to_list = False
        for m in range(M):
            if neighbour_array[m, atom] == None and not added_to_list:
                neighbour_array[m, atom] = nearest_sites_list[index]
                # print(f"added{nearest_sites_list[index]} to the list of neigbours of site{atom}")
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
    :param sampling_freq: Frequency (nb of steps) at which the position of the atoms should be sampled
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


def initialize_index_vector(structure: pmg.Structure) -> np.arry:
    """Initialize an index vector which links the occupancy vector(species at each site) and the current position of the atoms based on their starting site

    :param structure: supercell
    :type structure: pmg.Structure
    :return: 1D array with position of atoms at the initial situation
    :rtype: np.array
    """
    return np.arange(structure.num_sites)

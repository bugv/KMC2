"""This module contains the functions related to frequency calculations

It contains functions to initially attribute a frequency to hops with different types of atoms 
and to calculate the frequency of given events

Available functions:
- frequency_calculator"""

import numpy as np
import structure_management


def frquency_calculator(atom_key: dict) -> np.array:
    """Function that creats a dict that associates the frequency of a hop between a vacancy and an atom of each type

    :param atom_key: key relating each type of atom to a unique int (use atom_key_builder to create it)
    :type atom_key: dict
    :return: array where the value corresponds to the frequency of a hop for the atom with that index as identification
    :rtype: np.array
    """
    frequency_vector = np.full(len(atom_key), 1 / len(atom_key))
    return frequency_vector


def standardize_frequencies(user_frequencies: dict, atom_key: dict) -> np.array:
    """Function that takes the frequencies per atom inputed by the user and converts them to a freqeuency array

    :param user_frequencies: dict containing as keys the species and as values their frequency, the vacancy should be included as {"V0+":0}
    :type user_frequencies: dict
    :param atom_key: The atom key dict that relates the species in the structure with a unique int, generated with atom_key_builder
    :type atom_key: dict
    :return: 1D array with the frequency of a given atom at the index that corresponds to that atom in the array
    :rtype: np.array
    """
    frequency_vector = np.full(len(atom_key), None)
    for key in atom_key:
        frequency_vector[atom_key[key]] = user_frequencies[key]
    return frequency_vector


# NOTE The vacancy is the dummy species X, which is described as "X0+" in the structure


def hop_frequency_calculator(
    vacancy_position: int,
    occupancy_vector: np.array,
    neighbour_list: np.array,
    frequency_vect: np.array,
) -> tuple:
    """Function that calculates the frequency for the possible
       hops with the current position of the atoms

    :param vacancy_position: index indicating position of the vacancy in the structure
    :type vacancy_position: int
    :param occupancy_vector: occupancy vector
    :type occupancy_vector: np.array
    :param neighbour_list: list of neighbours
    :type neighbour_list: np.array
    :param frequency_vect: frequency of hops depending on atom type
    :type frequency_vect: np.array
    :return: frequency of the hops for the atoms that are neighbouring the vacancy
             and the sum of the frequencies
    :rtype: tuple
    """
    neighbours = neighbour_list[:, vacancy_position]
    neighbours = np.array([x for x in neighbours if x is not None])
    freq_neighbours = np.full(len(neighbours), None)
    for count, values in enumerate(neighbours):
        freq_neighbours[count] = frequency_vect[int(occupancy_vector[int(values)])]
    sum_frequencies = np.sum(freq_neighbours)
    freq_neighbours = (np.cumsum(freq_neighbours)) / sum_frequencies
    return freq_neighbours, sum_frequencies


def select_event(freq_neighbours: np.array) -> int:
    """Function to select an event from the frequency vector,
       returns the index of the neighbour with which the vacancy will switch

    :param freq_neighbours: 1D array with cumulative frequencies which range from 0 to 1
    :type freq_neighbours: np.array
    :return: an int corresponding to the column index of the neighbour the vacancy will switch with
    :rtype: int
    """
    rho = random_number()
    return np.min(np.where(freq_neighbours > rho))


def random_number() -> float:
    """Random number generator

    :return: Random number between 0 and 1
    :rtype: float
    """
    rng = np.random.default_rng()
    return rng.random()


def time_step_calculator(frequency_sum: float) -> float:
    """Function to calculate the time step

    :param frequency_sum: The sum of the freqeuency of all possible events
    :type frequency_sum: float
    :return: The time step for that event
    :rtype: float
    """
    return (-1) / (frequency_sum) * np.log(random_number())

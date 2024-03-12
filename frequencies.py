"""This module contains the functions related to frequency calculations

It contains functions to initially attribute a frequency to hops with different types of atoms 
and to calculate the frequency of given events

Available functions:
- frequency_calculator"""

import numpy as np


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

    :param user_frequencies: dict containing as keys the species and as values their frequency
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


## Hop frequency calculator function
# Input: occupancy vector, list of neighbours
# Output: vector with the event frequencies  (same order as the order of the neighbours in the list)
# 1. Find position of vacancy in the occupancy vector
# 2. get corresponding neighbours from the list of neighbours
# 3. create vector of size number of neighbours to stop the frequencies
# 4. For each neighbour in the list check with key the type of atom, then get frequency from list of frequencies as a function of atom type and put it into the vector
# 5. Calculate sum of all values in the frequency vector
# 6. divide frequency vector by sum
# 7. add previous value in vector to the value in vector (get end of interval in the sketch)
# 8. Return vector and sum


# def hop_frequency_calculator(occupancy_vector: np.array, neighbour_list:np.array)-> tuple:
# return


## Choose event function
# Input: random number rho, event frequency vector
# Output: event (number of the neighbour which is swapped?)
# 1. loop check if random number is greater than nth value in vector if no return n if yes check n+1th value...
# 2. return n (check if this makes sense)


def random_number() -> int:
    """Random number generator

    :return: Random number between 0 and 1
    :rtype: int
    """
    rng = np.random.default_rng()
    return rng.random()

"""This module contains the functions to execute the main parts of the KMC algorithm

It conatins functions corresponding to the initialization, 
the main body of the algorithm and the outputting of results

Available functions:
-Initialization"""

import structure_management
import frequencies
import calculators
import numpy as np
from structure_management import AtomPositions


def intialization(
    structure,
    user_frequencies: dict,
    radius: float,
    sampling_frequency: float,
    total_nb_steps: int,
) -> tuple:
    """Function that runs all of the steps of the initialization

    :param structure: supercell
    :type structure: _type_
    :param user_frequencies: dict, where the keys are the species in the structure and the values the associated frequency
    :type user_frequencies: dict
    :param radius: Cutoff radius for nearest neighbour
    :type radius: float
    :param sampling_frequency: Frequency (nb of steps) at which the position of the atoms should be sampled
    :type sampling_frequency: float
    :param total_nb_steps: total number of steps to perform
    :type total_nb_steps: int
    :return: Tuple containing in order the occupancy vector, the neighbour vector,
    the atom key, the frequency vector, the index_array,
    the current position array, the data collector array,
    the total number of steps, the sampling frequency and
    the current position of the vacancy
    :rtype: tuple
    """
    atom_key = structure_management.atom_key_builder(structure)
    occupancy_vector = structure_management.occupancy_vector_builder(
        structure, atom_key
    )
    neighour_array = structure_management.neighbour_finder(structure, radius)
    frequency_vector = frequencies.standardize_frequencies(user_frequencies, atom_key)
    current_position_array = structure_management.current_position_array_builder(
        structure
    )
    data_collector = structure_management.data_collector_builder(
        structure, sampling_frequency, total_nb_steps
    )
    index_array = structure_management.index_vector_builder(structure)
    vac_position = structure_management.find_vac(occupancy_vector, atom_key)
    atom_pos = AtomPositions(occupancy_vector, current_position_array, index_array)
    return (
        atom_pos,  # 0
        atom_key,  # 1
        neighour_array,  # 2
        frequency_vector,  # 3
        data_collector,  # 4
        total_nb_steps,  # 5
        sampling_frequency,  # 6
        vac_position,  # 7
    )


##Driver function
# This function implements the loop in the algorithm, it takes as input functions corresponding to different steps
# Input:  tuple from the initialize function, hop_frequency_calculator function, random_number_generator function, acceptor function
# Output:


def driver_v2(
    initialized_values: tuple,
) -> np.array:
    """Function that executes the main loop of the KMC

    :param initialized_values: Tuple returned by the initialization function
    :type initialized_values: tuple
    :return: data collector array
    :rtype: np.array
    """
    vac_position = initialized_values[7]
    t = 0
    r = initialized_values[4]
    atom_pos = initialized_values[0]
    for nb_steps in range(0, initialized_values[5]):
        freq_neighbours, sum_freq = frequencies.hop_frequency_calculator(
            vac_position,
            atom_pos.occupancy_vector,
            initialized_values[2],
            initialized_values[3],
            initialized_values[1],
        )
        event = frequencies.select_event(freq_neighbours)
        event = initialized_values[2][event, vac_position]
        t = t + frequencies.time_step_calculator(sum_freq)
        atom_pos.swap(vac_position, event)
        vac_position = event
        if nb_steps % initialized_values[6] == 0:
            r[:, :, int(nb_steps / initialized_values[6])] = (
                atom_pos.current_position_array
            )
    return r

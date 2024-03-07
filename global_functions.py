"""This module contains the functions to execute the main parts of the KMC algorithm

It conatins functions corresponding to the initialization, 
the main body of the algorithm and the outputting of results

Available functions:
-Initialization"""

import structure_management
import frequencies
import calculators


def intialization(
    structure,
    radius: float,
    sampling_frequency: float,
    total_nb_steps: int,
) -> tuple:
    """Function that runs all of the steps of the initialization

    :param structure: supercell
    :type structure: _type_
    :param radius: Cutoff radius for nearest neighbour
    :type radius: float
    :param sampling_frequency: Frequency (nb of steps) at which the position of the atoms should be sampled
    :type sampling_frequency: float
    :param total_nb_steps: total number of steps to perform
    :type total_nb_steps: int
    :return: Tuple containing the occupancy vector, the neighbour vector, the frequency vector,
    the current position array, the data collector array, the current time(0), and the number of steps performed(0)
    :rtype: tuple
    """
    atom_key = structure_management.atom_key_builder(structure)
    occupancy_vector = structure_management.occupancy_vector_builder(
        structure, atom_key
    )
    neighour_vector = structure_management.neighbour_finder(structure, radius)
    frequency_vector = frequencies.frquency_calculator(atom_key)
    current_position_array = structure_management.current_position_array_builder(
        structure
    )
    data_collector = structure_management.data_collector_builder(
        structure, sampling_frequency, total_nb_steps
    )
    t = 0
    nb_steps = 0
    return (
        occupancy_vector,
        neighour_vector,
        frequency_vector,
        current_position_array,
        data_collector,
        t,
        nb_steps,
    )


##Driver function
# This function implements the loop in the algorithm, it takes as input functions corresponding to different steps


##Result analysis function
# Thic function

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
    equivalent_sites_array = structure_management.get_equiv_in_primative(structure)
    displacements_tensor = structure_management.get_displacement_tensor(
        structure_management.empty_structure(structure), radius
    )
    neighour_array = structure_management.get_neighbours_from_displ(
        structure_management.empty_structure(structure),
        equivalent_sites_array,
        radius,
        displacements_tensor,
    )
    frequency_vector = frequencies.standardize_frequencies(user_frequencies, atom_key)
    current_position_array = structure_management.current_position_array_builder(
        structure
    )
    data_collector = structure_management.data_collector_builder(
        structure, sampling_frequency, total_nb_steps
    )
    index_array = structure_management.index_vector_builder(structure)
    vac_position = structure_management.find_vac(occupancy_vector, atom_key)
    atom_pos = AtomPositions(
        occupancy_vector,
        current_position_array,
        index_array,
        neighour_array,
        equivalent_sites_array,
        displacements_tensor,
    )
    time_collector = structure_management.time_collector_builder(
        sampling_frequency, total_nb_steps
    )
    return {
        "atom_pos": atom_pos,
        "atom_key": atom_key,
        "equivalent_sites_array": equivalent_sites_array,
        "neighbour_array": neighour_array,
        "frequency_vector": frequency_vector,
        "data_collector": data_collector,
        "total_nb_steps": total_nb_steps,
        "sampling_frequency": sampling_frequency,
        "vac_position": vac_position,
        "time_collector": time_collector,
        "time": 0,
        "displacement_tensor": displacements_tensor,
    }
    # return (
    #     atom_pos,  # 0
    #     atom_key,  # 1
    #     neighour_array,  # 2
    #     frequency_vector,  # 3
    #     data_collector,  # 4
    #     total_nb_steps,  # 5
    #     sampling_frequency,  # 6
    #     vac_position,  # 7
    # )


##Driver function
# This function implements the loop in the algorithm, it takes as input functions corresponding to different steps
# Input:  tuple from the initialize function, hop_frequency_calculator function, random_number_generator function, acceptor function
# Output:


def driver(
    initialized_values: dict,
) -> np.array:
    """Function that executes the main loop of the KMC

    :param initialized_values: Tuple returned by the initialization function
    :type initialized_values: tuple
    :return: data collector array
    :rtype: np.array
    """
    vac_position = initialized_values["vac_position"]
    data_collector = initialized_values["data_collector"]
    atom_pos = initialized_values["atom_pos"]
    for nb_steps in range(0, initialized_values["total_nb_steps"]):
        freq_neighbours, sum_freq = frequencies.hop_frequency_calculator(
            initialized_values["vac_position"],
            initialized_values["atom_pos"].occupancy_vector,
            initialized_values["neighbour_array"],
            initialized_values["frequency_vector"],
            initialized_values["atom_key"],
        )
        event = frequencies.select_event(freq_neighbours)
        event = initialized_values["neighbour_array"][
            event, initialized_values["vac_position"]
        ]
        initialized_values["time"] = initialized_values[
            "time"
        ] + frequencies.time_step_calculator(sum_freq)
        initialized_values["atom_pos"].swap(initialized_values["vac_position"], event)
        initialized_values["vac_position"] = event
        if nb_steps % initialized_values["sampling_frequency"] == 0:
            initialized_values["data_collector"][
                :, :, int(nb_steps / initialized_values["sampling_frequency"])
            ] = initialized_values["atom_pos"].current_position_array
    return initialized_values

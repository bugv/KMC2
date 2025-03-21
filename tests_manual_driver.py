import structure_management
import frequencies
import numpy as np
from pymatgen.core import Structure
import pymatgen.core as pmg
from structure_management import AtomPositions
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar
import json
import copy
import global_functions


"""Function that executes the main loop of the KMC

:param initialized_values: Tuple returned by the initialization function
:type initialized_values: tuple
:return: data collector array
:rtype: np.array
"""

initialized_values = global_functions.read_full_from_json("Initialized_struct.json")
nb_steps = initialized_values["total_nb_steps"]

random_array = np.random.random((2, initialized_values["total_nb_steps"]))
random_array = np.zeros((2, initialized_values["total_nb_steps"]))

equal_freqs = True
# if yes -> calculate frequencies once, set a variable to true (fast)
# if no -> set the variale to false (slower)
if equal_freqs:
    freq_neighbours, sum_freq = frequencies.hop_frequency_calculator(
        initialized_values["vac_position"],
        initialized_values["atom_pos"].occupancy_vector,
        initialized_values["neighbour_array"],
        initialized_values["frequency_vector"],
    )

def single_step(nb_steps,initialized_values) :
    nb_steps = 0
    event = frequencies.select_event_alternative(
        random_array, nb_steps, freq_neighbours
    )
    print(event)

    initialized_values["time"] = initialized_values[
        "time"
    ] + frequencies.time_step_calculator_alternative(
        random_array, nb_steps, sum_freq
    )

    initialized_values["atom_pos"].swap_displ_V2(
        initial=initialized_values["vac_position"], position_final_in_list=event
    )
    event = initialized_values["neighbour_array"][
        event, initialized_values["vac_position"]
    ]
    print(event)

    initialized_values["vac_position"] = event
    if nb_steps % initialized_values["sampling_frequency"] == 0:

        initialized_values["data_collector"][
            :, :, int(nb_steps / initialized_values["sampling_frequency"])
        ] = initialized_values["atom_pos"].current_position_array
        initialized_values["time_collector"][
            0, int(nb_steps / initialized_values["sampling_frequency"])
        ] = initialized_values["time"]



single_step(0,initialized_values); 
print(initialized_values["vac_position"])
print(np.where(initialized_values["atom_pos"].occupancy_vector == 2))
initialized_values["atom_pos"].current_position_array[:,initialized_values["vac_position"]]
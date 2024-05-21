"""This module contains the functions to execute the main parts of the KMC algorithm

It conatins functions corresponding to the initialization, 
the main body of the algorithm and the outputting of results

Available functions:
-Initialization
-Driver"""

import structure_management
import frequencies
import calculators
import numpy as np
from pymatgen.core import Lattice, Structure
import pymatgen.core as pmg
from structure_management import AtomPositions
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import json
import time
import copy
import psutil


def read_input_file(file_name: str) -> tuple:
    """Reads the data from an input json file and returns a tuple with all of the data.
       Data is ordered such that it can be directly input into the initialization function

    :param file_name: path to the input file (json)
    :type file_name: str
    :return: input ready for initialization, in order: supercell, frequency dict,
             neighbour radius, sampling frequency, number of steps
    :rtype: tuple
    """
    with open(file_name) as json_file:
        data = json.load(json_file)
        number_steps = data["number_steps"]
        samping_frequency = data["sampling_frequency"]
        neighbour_radius = data["neighbour_radius"]
        frequencies_dict = data["hop_frequencies"]
        frequencies_dict["X0+"] = 0
        lattice = data["lattice"]
        supercell_size = data["supercell_size"]
        atoms = data["atoms"]
        coords = data["coords"]
        struct = Structure(lattice, atoms, coords)
        primcell = SpacegroupAnalyzer(
            structure_management.empty_structure(struct)
        ).find_primitive()
        struct = Structure(lattice, atoms, coords)
        supercell = struct * supercell_size
        return (
            primcell,  # 0
            supercell,  # 1
            frequencies_dict,  # 2
            neighbour_radius,  # 3
            samping_frequency,  # 4
            number_steps,  # 5
        )


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def write_full_to_json(full_struct: dict, file_name: str):
    full_struct = copy.deepcopy(full_struct)
    full_struct["occupancy_vector"] = full_struct["atom_pos"].occupancy_vector
    full_struct["index_array"] = full_struct["atom_pos"].index_array
    full_struct["neighbour_array"] = full_struct["atom_pos"].neighbour_list
    full_struct["equivalent_sites_in_prim"] = full_struct[
        "atom_pos"
    ].equivalent_sites_in_prim.tolist()
    full_struct["current_position_array"] = full_struct[
        "atom_pos"
    ].current_position_array.tolist()
    full_struct["displacement_tensor"] = full_struct["atom_pos"].displacement_tensor
    del full_struct["atom_pos"]
    with open(file_name, "w") as json_file:
        json.dump(full_struct, json_file, cls=NumpyEncoder)


def read_full_from_json(file_name: str) -> dict:
    with open(file_name, "r") as json_file:
        full_struct = json.load(json_file)
        full_struct["atom_pos"] = AtomPositions(
            np.asarray(full_struct["occupancy_vector"]),
            np.asarray(full_struct["current_position_array"]),
            np.asarray(full_struct["index_array"]),
            np.asarray(full_struct["neighbour_array"]),
            np.asarray(full_struct["equivalent_sites_in_prim"]),
            np.asarray(full_struct["displacement_tensor"]),
        )
        del full_struct["occupancy_vector"]
        del full_struct["current_position_array"]
        del full_struct["index_array"]
        del full_struct["equivalent_sites_in_prim"]
        full_struct["equivalent_sites_array"] = np.asarray(
            full_struct["equivalent_sites_array"]
        )
        full_struct["neighbour_array"] = np.asarray(full_struct["neighbour_array"])
        full_struct["data_collector"] = np.asarray(full_struct["data_collector"])
        full_struct["frequency_vector"] = np.asarray(full_struct["frequency_vector"])
        full_struct["time_collector"] = np.asarray(full_struct["time_collector"])
        full_struct["displacement_tensor"] = np.asarray(
            full_struct["displacement_tensor"]
        )

    return full_struct


def initialization(
    primcell,
    supercell,
    user_frequencies: dict,
    radius: float,
    sampling_frequency: float,
    total_nb_steps: int,
) -> tuple:
    """Function that runs all of the steps of the initialization

    :param primcell: primitive cell
    :type primcell: pmg.Structure
    :param supercell: supercell
    :type structure: pmg.Structure
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
    supercell, initial_vac_pos = structure_management.add_vacancy_random(supercell)
    # supercell = structure_management.add_vacancy(supercell, initial_vac_pos)

    start_time = time.time()
    process = psutil.Process()
    mem_before = process.memory_info().rss / 1024  # in kilobytes
    atom_key = structure_management.atom_key_builder(supercell)
    mem_after = process.memory_info().rss / 1024  # in kilobytes
    mem_used = mem_after - mem_before
    print("mem used for atom key builder (kilobytes)", mem_used)
    end_time = time.time()
    print("Atom key builder time", end_time - start_time)

    start_time = time.time()
    process = psutil.Process()
    mem_before = process.memory_info().rss / 1024  # in kilobytes
    occupancy_vector = structure_management.occupancy_vector_builder(
        supercell, atom_key
    )
    mem_after = process.memory_info().rss / 1024  # in kilobytes
    mem_used = mem_after - mem_before
    print("mem used for occupancy vector builder(kilobytes)", mem_used)
    end_time = time.time()
    print("Occupancy vector builder time", end_time - start_time)

    start_time = time.time()
    process = psutil.Process()
    mem_before = process.memory_info().rss / 1024  # in kilobytes
    equivalent_sites_array = structure_management.get_equiv_in_primative(
        primcell, supercell
    )
    mem_after = process.memory_info().rss / 1024  # in kilobytes
    mem_used = mem_after - mem_before
    print("mem used for get equiv in primative (kilobytes)", mem_used)
    end_time = time.time()
    print("Get equivalent in prim time", end_time - start_time)

    start_time = time.time()
    process = psutil.Process()
    mem_before = process.memory_info().rss / 1024  # in kilobytes
    displacements_tensor = structure_management.get_displacement_tensor(
        primcell, supercell, radius
    )
    mem_after = process.memory_info().rss / 1024  # in kilobytes
    mem_used = mem_after - mem_before
    print("mem used for get displacements (kilobytes)", mem_used)
    end_time = time.time()
    print("get displacement tensor time", end_time - start_time)

    start_time = time.time()
    process = psutil.Process()
    mem_before = process.memory_info().rss / 1024  # in kilobytes
    neighour_array = structure_management.alternate_neighbour_finder(
        supercell,
        radius,
    )
    mem_after = process.memory_info().rss / 1024  # in kilobytes
    mem_used = mem_after - mem_before
    print("mem used for get neighbours (kilobytes)", mem_used)
    end_time = time.time()
    print("Get neighbour array time", end_time - start_time)

    start_time = time.time()
    process = psutil.Process()
    mem_before = process.memory_info().rss / 1024  # in kilobytes
    frequency_vector = frequencies.standardize_frequencies(user_frequencies, atom_key)
    current_position_array = structure_management.current_position_array_builder(
        supercell
    )
    mem_after = process.memory_info().rss / 1024  # in kilobytes
    mem_used = mem_after - mem_before
    print("mem used for get freq vect (kilobytes)", mem_used)
    end_time = time.time()
    print("Get frequency vect time", end_time - start_time)

    start_time = time.time()
    process = psutil.Process()
    mem_before = process.memory_info().rss / 1024  # in kilobytes
    data_collector = structure_management.data_collector_builder(
        supercell, sampling_frequency, total_nb_steps
    )
    mem_after = process.memory_info().rss / 1024  # in kilobytes
    mem_used = mem_after - mem_before
    print("mem used for build data collector (kilobytes)", mem_used)
    end_time = time.time()
    print("build data collector time", end_time - start_time)

    mem_before = process.memory_info().rss / 1024  # in kilobytes
    index_array = structure_management.index_vector_builder(supercell)

    end_time = time.time()
    print("index vector builder time", end_time - start_time)

    start_time = time.time()
    process = psutil.Process()
    mem_before = process.memory_info().rss / 1024  # in kilobytes
    vac_position = structure_management.find_vac(occupancy_vector, atom_key)
    mem_after = process.memory_info().rss / 1024  # in kilobytes
    mem_used = mem_after - mem_before
    print("mem used for find vac (kilobytes)", mem_used)
    end_time = time.time()
    print("get vac pos time", end_time - start_time)

    atom_pos = AtomPositions(
        occupancy_vector,
        current_position_array,
        index_array,
        neighour_array,
        equivalent_sites_array,
        displacements_tensor,
    )

    start_time = time.time()
    time_collector = structure_management.time_collector_builder(
        sampling_frequency, total_nb_steps
    )
    end_time = time.time()
    print("time collector builder time", end_time - start_time)
    print("Completed Initialization")
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
        "initial_vac_position": vac_position,
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
    # loop for number of steps required
    for nb_steps in range(0, initialized_values["total_nb_steps"]):
        # get the frequencies of a swap with the neighbours of the vacancy
        freq_neighbours, sum_freq = frequencies.hop_frequency_calculator(
            initialized_values["vac_position"],
            initialized_values["atom_pos"].occupancy_vector,
            initialized_values["neighbour_array"],
            initialized_values["frequency_vector"],
        )
        # select the swap -> index of neighbour with which vacancy will switch
        event = frequencies.select_event(freq_neighbours)
        # print("selected event which neighbour", event)
        # select swap get index of swap site in structure
        event = initialized_values["neighbour_array"][
            event, initialized_values["vac_position"]
        ]
        # print("selected event actual index of neighbour", event)
        # update time
        initialized_values["time"] = initialized_values[
            "time"
        ] + frequencies.time_step_calculator(sum_freq)
        # swap
        # print(
        #     "current position before swap",
        #     initialized_values["atom_pos"].current_position_array,
        # )
        initialized_values["atom_pos"].swap_displ(
            initial=initialized_values["vac_position"], final=event
        )
        # print(
        #     "current position after swap",
        #     initialized_values["atom_pos"].current_position_array,
        # )
        # update vacancy position
        initialized_values["vac_position"] = event
        if nb_steps % initialized_values["sampling_frequency"] == 0:
            # print(
            #     "current position when sampled",
            #     initialized_values["atom_pos"].current_position_array,
            # )
            initialized_values["data_collector"][
                :, :, int(nb_steps / initialized_values["sampling_frequency"])
            ] = initialized_values["atom_pos"].current_position_array
            initialized_values["time_collector"][
                0, int(nb_steps / initialized_values["sampling_frequency"])
            ] = initialized_values["time"]
            # print("time recorded", initialized_values["time_collector"])
            # print(
            #     "when included in data collector",
            #     initialized_values["data_collector"][
            #         :, :, int(nb_steps / initialized_values["sampling_frequency"])
            #     ],
            # )
    print("Completed run")
    return initialized_values

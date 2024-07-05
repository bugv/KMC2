"""This module contains the functions to execute the main parts of the KMC algorithm

It conatins functions corresponding to the initialization, 
the main body of the algorithm and the outputting of results

Available functions:
-read_input_file
-read_input_file_composition
-write_full_to_json
-write_struct_to_poscar
-initialization
-driver"""

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
        frequencies_dict["X0+"] = 0.0
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
            primcell,
            supercell,
            frequencies_dict,
            neighbour_radius,
            samping_frequency,
            number_steps,
        )


def read_input_file_composition(file_name: str) -> tuple:
    """Reads the data from an input json file and returns a tuple with all of the data.
       Data is ordered such that it can be directly input into the initialization function
       This version is to crate a random structure given  a concentration

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
        frequencies_dict["X0+"] = 0.0
        lattice = data["lattice"]
        supercell_size = data["supercell_size"]
        atoms = data["atoms"]
        coords = data["coords"]
        composition = data["composition"]
        struct = Structure(lattice, atoms, coords)
        primcell = SpacegroupAnalyzer(
            structure_management.empty_structure(struct)
        ).find_primitive()
        struct = Structure(lattice, atoms, coords)
        supercell = struct * supercell_size
        supercell = structure_management.create_alloy(supercell, composition)
        return (
            primcell,
            supercell,
            frequencies_dict,
            neighbour_radius,
            samping_frequency,
            number_steps,
        )


class NumpyEncoder(json.JSONEncoder):
    """encoder used to manage writing numpy arrays to json files"""

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def write_full_to_json(full_struct: dict, file_name: str):
    """Function to output a dict containing a full structure (AtomPosition object  and other arrays) to a json file

    :param full_struct: contains AtomPosition object and other arrays
    :type full_struct: dict
    :param file_name: name output json file
    :type file_name: str
    """
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
    """Function to read a full structure from a json file and convert it back to the corresponding arrays/ AtomPosition object

    :param file_name: name of json file
    :type file_name: str
    :return: Dict containing a AtomPosition object and the other arrays
    :rtype: dict
    """
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


def write_struct_to_poscar(structure: pmg.Structure) -> None:
    """Function to write a pymatgen structure as a POSCAR file. The output will have the name "POSCAR"

    :param structure: the structure to write to the file
    :type structure: pmg.Structure
    """
    structure_as_poscar = Poscar(structure, sort_structure=True)
    structure_as_poscar.write_file("POSCAR")


def initialization(
    primcell,
    supercell,
    user_frequencies: dict,
    radius: float,
    sampling_frequency: float,
    total_nb_steps: int,
) -> tuple:
    """Function that runs all of the steps of the initialization and writes the initial supercell to a POSCAR file

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

    atom_key = structure_management.atom_key_builder(supercell)

    compositon_dict = structure_management.get_comp_dict(supercell, atom_key)

    occupancy_vector = structure_management.occupancy_vector_builder(
        supercell, atom_key
    )

    equivalent_sites_array = structure_management.get_equiv_in_primative(
        primcell, supercell
    )

    displacements_tensor = structure_management.get_displacement_tensor(
        primcell, supercell, radius
    )

    neighour_array = structure_management.alternate_neighbour_finder(
        supercell,
        radius,
    )

    frequency_vector = frequencies.standardize_frequencies(user_frequencies, atom_key)
    current_position_array = structure_management.current_position_array_builder(
        supercell
    )

    data_collector = structure_management.data_collector_builder(
        supercell, sampling_frequency, total_nb_steps
    )

    index_array = structure_management.index_vector_builder(supercell)

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

    write_struct_to_poscar(supercell)
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
        "initial_vac_position": initial_vac_pos,
        "composition_dict": compositon_dict,
    }


def driver(
    initialized_values: dict,
) -> np.array:
    """Function that executes the main loop of the KMC

    :param initialized_values: Tuple returned by the initialization function
    :type initialized_values: tuple
    :return: data collector array
    :rtype: np.array
    """

    random_array = np.random.random((2, initialized_values["total_nb_steps"]))

    ### Accelerate if the frequencies are the same
    frequency_vect = initialized_values["frequency_vector"].astype(float)
    equal_freqs = np.allclose(
        frequency_vect[:-1],
        frequency_vect[0],
        rtol=1e-5,
        atol=1e-8,
    )
    if equal_freqs:
        freq_neighbours, sum_freq = frequencies.hop_frequency_calculator(
            initialized_values["vac_position"],
            initialized_values["atom_pos"].occupancy_vector,
            initialized_values["neighbour_array"],
            initialized_values["frequency_vector"],
        )
    # if yes -> calculate frequencies once, set a variable to true (fast)
    # if no -> set the variale to false (slower)

    for nb_steps in range(0, initialized_values["total_nb_steps"]):

        if not equal_freqs:
            freq_neighbours, sum_freq = frequencies.hop_frequency_calculator(
                initialized_values["vac_position"],
                initialized_values["atom_pos"].occupancy_vector,
                initialized_values["neighbour_array"],
                initialized_values["frequency_vector"],
            )

        event = frequencies.select_event_alternative(
            random_array, nb_steps, freq_neighbours
        )

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

        initialized_values["vac_position"] = event
        if nb_steps % initialized_values["sampling_frequency"] == 0:

            initialized_values["data_collector"][
                :, :, int(nb_steps / initialized_values["sampling_frequency"])
            ] = initialized_values["atom_pos"].current_position_array
            initialized_values["time_collector"][
                0, int(nb_steps / initialized_values["sampling_frequency"])
            ] = initialized_values["time"]

    print("Completed run")
    return initialized_values

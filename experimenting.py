import structure_management
import frequencies
import global_functions
import pymatgen.core as pmg
from pymatgen.core import Lattice, Structure
import numpy as np

from structure_management import AtomPositions

bcc_fe_Al = Structure(Lattice.cubic(2.8), ["Fe", "Al"], [[0, 0, 0], [0.5, 0.5, 0.5]])
polonium = Structure(Lattice.cubic(3.4), ["Po"], [[0.0, 0.0, 0.0]])
test_supercell = bcc_fe_Al * (2, 2, 2)
test_supercell2 = polonium * (3, 3, 3)

# print(bcc_fe_Al)

# test_supercell[3] =

# test = comp.get_el_amt_dict()

# print(bcc_fe_Al.species) #get a list eg. [Element Fe, Element Al]

# print(bcc_fe_Al.species)

# print(type(test_supercell))

# print(test_supercell.as_dataframe())

# print(test_supercell.as_dict(verbosity=0, fmt="abivars"))
# print(test_supercell)
# print(type(test_supercell.species))

# print(test_supercell[0])
# print(test_supercell[1])
# print(test_supercell[2])
# print(test_supercell[3])
# print(test_supercell[4])
# print(test_supercell[30].species.elements)
# print(test_supercell.num_sites)


# print(type(test_supercell.elements))

# print(structure_management.atom_key_builder(test_supercell))

test_supercell[5] = pmg.DummySpecies("X")
atom_key = structure_management.atom_key_builder(test_supercell)
user_freq = {"Al": 0.25, "Fe": 0.75, "X0+": 0}

initialized_structure = global_functions.intialization(
    test_supercell, user_freq, 2.5, 10, 100
)

print(global_functions.driver(initialized_structure))


# print(initialized_structure)
# print("occupancy vector", initialized_structure[0])

# print(initialized_structure[2])

# print(initialized_structure[1])


# print(structure_management.find_vac(initialized_structure[0], initialized_structure[1]))

# print(initialized_structure[1])
# frequencies.hop_frequency_calculator(
#     initialized_structure[0],
#     initialized_structure[2],
#     initialized_structure[3],
#     initialized_structure[1],
# )

# print(
#     frequencies.hop_frequency_calculator(
#         initialized_structure[0],
#         initialized_structure[2],
#         initialized_structure[3],
#         initialized_structure[1],
#     )
# )

# pmg.DummySpecies("X")


# test_atom_pos = AtomPositions(
#     initialized_structure[0], initialized_structure[5], initialized_structure[4]
# )
# print(test_atom_pos.current_position_array)
# test_atom_pos.swap(5, 6)
# print(test_atom_pos.current_position_array)


# print(
#     frequencies.hop_frequency_calculator(
#         initialized_structure[0],
#         initialized_structure[2],
#         initialized_structure[3],
#         initialized_structure[1],
#     )
# )

# freqtest = frequencies.hop_frequency_calculator(
#     initialized_structure[0],
#     initialized_structure[2],
#     initialized_structure[3],
#     initialized_structure[1],
# )[0]
# print(freqtest)
# rho = frequencies.random_number()
# print(rho)
# rho = 0.9
# print(np.where(freqtest > rho))
# print(np.min(np.where(freqtest > rho)))

# print(len(atom_key))
# print(frequencies.frquency_calculator(atom_key))
# print(frequencies.random_number())
# print(frequencies.random_number())


# structure_management.neighbour_finder_v2(test_supercell, 2.5)


# m = 8
# n = test_supercell.num_sites
# neighbour_list = test_supercell.get_neighbor_list(2.5)
# centers_list = neighbour_list[0]
# nearest_sites_list = neighbour_list[1]
# unique, counts = np.unique(centers_list, return_counts=True)
# counts[2] = 5
# print(counts)


# position_array = np.full((m, n), None, dtype=object)
# position_array = np.tile(np.arange(m), (n, 1)).T

# print(position_array)

# values = np.arange(m)
# position_array = np.where(position_array[:, :] > counts[:], values, None)


# if condition true use values if false use none
# test_initialize = global_functions.intialization(test_supercell, 2.5, 5, 100)
# print(type(test_initialize), len(test_initialize))

# print(test_initialize[1])
# print(structure_management.occupancy_vector_builder(test_supercell, atom_key))

import structure_management
import pymatgen.core as pmg

import numpy as np

structure = pmg.Structure(pmg.Lattice.cubic(2.8), ["Fe", "Al"],[[0,0,0], [0.5,0.5,0.5]]) * (2,2,2)

radius  =2.5

def test_fill_neighbour_array():
    neighbour_list = structure.get_neighbor_list(radius)[1]
    neighbour_array = structure_management.neighbour_finder(structure, radius)
    assert np.array_equal(neighbour_list, neighbour_array.flatten('F'))



import pytest
import pymatgen.core as pmg
import numpy as np
import structure_management


@pytest.fixture
def structure():
    return pmg.Structure(
        pmg.Lattice.cubic(2.8), ["Fe", "Al"], [[0, 0, 0], [0.5, 0.5, 0.5]]
    ) * (2, 2, 2)


@pytest.fixture
def radius():
    return 2.5


def test_fill_neighbour_array(structure, radius):
    neighbour_list = structure.get_neighbor_list(radius)[1]
    neighbour_array = structure_management.neighbour_finder(structure, radius)
    assert np.array_equal(neighbour_list, neighbour_array.flatten("F"))

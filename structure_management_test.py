import pytest
import pymatgen.core as pmg
import numpy as np
import structure_management


# Test get neighbour array (first version)
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


# Test max nb neighbours
@pytest.fixture
def struct_radius():
    bcc_struct = pmg.Structure(
        pmg.Lattice.cubic(2.8), ["Fe", "Al"], [[0, 0, 0], [0.5, 0.5, 0.5]]
    ) * (2, 2, 2)
    omega_struct = pmg.Structure(
        [
            [0.00000000, 0.00000000, 3.14710742],
            [2.52088978, 4.36630680, 0.00000000],
            [-5.04177956, 0.00000000, 0.00000000],
        ],
        ["Ti", "Ti", "Ti"],
        [
            [0, 0, 0],
            [0.5, 0.33333300, -0.33333300],
            [0.50000000, 0.66666600, -0.66666600],
        ],
    ) * (2, 2, 2)
    return [((bcc_struct, 2.5), 8), ((omega_struct, 4), 14)]


def test_max_nb_neighbours(struct_radius):
    for input, expected_r in struct_radius:
        structure, radius = input
        assert (
            structure_management.get_max_nb_neighbours(structure, radius) == expected_r
        )

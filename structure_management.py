## This module contains all the functions that mangage the structure
## and keeping track of the atoms

import numpy as np
import pymatgen.core as pmg
from pymatgen.core import Structure, Lattice

## Read from file to structure
# FUnction to read structure from file and produce a pymatgen structure
#Input: a file location
#Output: a pymatgen structure

## Create the occupancy vector
# Input: structure
# Output: Occupation vector of the structure with one vacancy in the structure (size=nb of sites in the structure)
# 1. read structure from the supercell
# 4. Get number of sites in the structure
# 5. Create a vector with length number of sites in the structure
# 6. Convert structure to occupancy vector by reading each item in the structure, and using key to get appropriate value to write to the vector
# 7. Return occupancy vector

def occupancy_vector_builder (structure: pmg.Structure, atom_key: dict) -> np.array:
    N = structure.num_sites # get number of sites in structure
    occupancy_vector = np.empty(N) #create empty occupancy vector of length  = nb of sites
    for i in range(N): #can't enumerate over empty array
        occupancy_vector[i] = atom_key[str(structure[i].species.elements[0])] #index 0 because strcuture[i].species.elements is a list with only one entry
    return occupancy_vector
    


##Create key that links each type of atom in the structure to an integer
#input: the structure
#output: a dict where the keys are the type of element and the values a corresponding integer eg: Vacancy->0, atom A->1, atom B->2 


def atom_key_builder (structure: pmg.Structure) -> dict:
    list_values = list(range(len(structure.elements)))
    element_list = structure.elements
    for i in range(len(element_list)): #NOTE: not ideal to use range of len but I couldn't enumerate over elements in list to convert them to string
        element_list[i] = str(element_list[i])
    atom_key = dict(zip(element_list, list_values))
    return atom_key
# TODO check if possible to use species not element -> might be better to include vacancies



##Creates a list of neighbours for each site
#Input: pymatgen structure read from the file 
#Output: Neighbour array(size = nb of sites x nb of neighbours per site(dependent on type of structure))
# 1. Call pymatgen get_all_neighbours need cut of distance for neighbours
# 2. Get number of neighbours per site to initialize array properly
# 2. Convert result to array
# 3. Return array

def neighbour_finder (structure: pmg.Structure, radius: float) -> np.array:
    N = structure.num_sites # Get number of sites N
    neighbour_list = structure.get_neighbor_list(radius) #NOTE cutoff value based on radius (but the value of radius is arbitrary, I used 2.5)
    centers_list = neighbour_list[0]
    nearest_sites_list = neighbour_list[1]
    print("nearest_sites_list", nearest_sites_list)
    #Get number of neighbours per site based on the number of neighbours of the first site
    M = 0
    for index, atom in enumerate(centers_list):
        if centers_list[index] ==0:
            M+=1
    # Create neighbours array and fill it
    neighbour_array = np.full((M, N), None)    
    for index, atom in enumerate(centers_list):
        added_to_list = False
        for m in range(M):
            if neighbour_array[m, atom] == None and not added_to_list:
                neighbour_array[m, atom] = nearest_sites_list[index]
                # print(f"added{nearest_sites_list[index]} to the list of neigbours of site{atom}")
                added_to_list = True
    return neighbour_array


##Initialize a temporary position array to store 
#Input: pymatgen stucture
#Output: array with size 3 x nb of sites (same structure as a single time slice of the R vector)
# 1. Use the get cartesian coordinates function from pymatgen to get the coordinates of each site and create array from them
# 2. return the array 

def temp_position_array_builder (structure: pmg.structure) -> np.array:
    return temp_position_array

## Initialize R, the data collection array
#Input: temorary position array,  number of times sampled/ sampling frequency? (include default value for case when it is not specified),
#Output:  array R (size: 3 x nb of sites x nb of  times sampled), contains initial positions of atoms, rest is empty
# 1. create array of size 3 x nb of sites x nb of times sampled
# 2. fill first time step layer with content of the temporary position array


def data_collector_builder (temp_position_array: np.array, sampling_freq: float) -> np.array:
    return data_collector


## Update temporary position array
#Input: event, temporary position array, list of neighbours, (position of vacancy ??)
#Output:  temporary position array

## Update structure (occupancy vector, temp position array and index vector)
#Input:
#Output:


## Initialize index vector (links the occupancy vector and the temp position array)
#Input:
#Output:


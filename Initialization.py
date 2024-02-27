## This module contains the functions that are run in initialization process
import numpy as np


## Create the occupancy vector
# Input: structure, atom_key
# Output: Occupation vector of the structure with one vacancy in the structure (size=nb of sites in the structure)
# 1. read structure from the supercell
# 2. Check number of vacancies in the supercell
#       a. 1 vacancy -> fine to continue
#       b. 0 vacancies -> Replace random atom in supercell with a vacancy
#       c. more than one vacancy -> Raise error? or replace all vacancies but one with atoms
# 4. Get number of sites in the structure
# 5. Create a vector with length number of sites in the structure
# 6. Convert structure to occupancy vector by reading each item in the structure, and using key to get appropriate value to write to the vector
# 7. Return occupancy vector

def occupancy_vector_builder (file_name: str) -> np.array:
    return occupancy_vector


##Create key that links each type of atom in the structure to an integer
#input: the structure
#output: a key that associates each element in the structure with an integer eg: Vacancy->0, atom A->1, atom B->2 (as an array? dict?)
#
def atom_key_builder (structure: pymatgen.structure) -> dict:
    return atom_key



##Creates a list of neighbours for each site
#Input: pymatgen structure read from the file 
#Output: Neighbour array(size = nb of sites x nb of neighbours per site(dependent on type of structure))
# 1. Call pymatgen get_all_neighbours need cut of distance for neighbours
# 2. Convert result to array
# 3. Return array

def neighbour_finder (structure: pymatgen.structure) -> np.array:
    return neighbour_array



##Frequency Calculator2. Check number of vacancies in the supercell
#       a. 1 vacancy -> fine to continue
#       b. 0 vacancies -> Replace random atom in supercell with a vacancy
#       c. more than one vacancy -> Raise error? or replace all vacancies but one with atoms
#Input: Key that relates the different types of atoms to ints
#Output: Vector with frequencies for each element in the structure (size: nb of different elements in the structure)
# 1. Create vector with length equal to number of elements in the structure 
# 2. return vector

def frquency_calculator (atom_key: dict)-> np.array:
    return frequency_vector


##Initialize a  position array to store position at current time t
#Input: pymatgen stucture
#Output: array with size 3 x nb of sites (same structure as a single time slice of the R vector)
# 1. Use the get cartesian coordinates function from pymatgen to get the coordinates of each site and create array from them
# 2. return the array 

def position_array_builder (structure: pymatgen.structure) -> np.array:
    return position_array


## Initialize R, the data collection array
#Input: temorary position array,  number of times sampled/ sampling frequency? (include default value for case when it is not specified),
#Output:  array R (size: 3 x nb of sites x nb of  times sampled), contains initial positions of atoms, rest is empty
# 1. create array of size 3 x nb of sites x nb of times sampled
# 2. fill first time step layer with content of the temporary position array


def data_collector_builder (temp_position_array: np.array, sampling_freq: float) -> np.array:
    return data_collector


## Initialization function that calls the other functions and sets appropriate variables to 0  
# -> should this be a function or should it just be included as such in the main file
#input: poscar file with super cell, sampling frequency
#output: array R, occupation vector, neighbour vector, frequency vector, nbofsteps=0, time=0
# 1. Read structre from poscar file and create a pymatgen structure
# 2. Call the key creation function ->returns a key
# 3. Call the create occupancy vector function-> return the occupancy vector
# 4. Call the create neighbour list function -> return the neighbour list array
# 5. Call the frequency calculator function -> returns a vector with the frequency for each element in the structure
# 6. Call the initialize temp position array -> returns an array with initial position
# 7. Call the initialize R function -> returns datacollection array R
# 8. Set time = 0
# 9. Set number of steps performed = 0


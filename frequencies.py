##Frequency Calculator
#Input: Key that relates the different types of atoms to ints
#Output: Vector with frequencies for each element in the structure (size: nb of different elements in the structure)
# 1. Create vector with length equal to number of elements in the structure 
# 2. return vector

def frquency_calculator (atom_key: dict)-> np.array:
    return frequency_vector


## Hop frequency calculator function
#Input: occupancy vector, list of neighbours
#Output: vector with the event frequencies  (same order as the order of the neighbours in the list)
# 1. Find position of vacancy in the occupancy vector
# 2. get corresponding neighbours from the list of neighbours
# 3. create vector of size number of neighbours to stop the frequencies
# 4. For each neighbour in the list check with key the type of atom, then get frequency from list of frequencies as a function of atom type and put it into the vector
# 5. Calculate sum of all values in the frequency vector
# 6. divide frequency vector by sum 
# 7. add previous value in vector to the value in vector (get end of interval in the sketch)
# 8. Return vector and sum


## Choose event function
#Input: random number rho, event frequency vector
#Output: event (number of the neighbour which is swapped?)
# 1. loop check if random number is greater than nth value in vector if no return n if yes check n+1th value...
# 2. return n (check if this makes sense)
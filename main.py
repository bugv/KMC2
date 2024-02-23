##Get as input from user:
#       - poscar file containing the supercell
#       - sampling frequency (optional, default value otherwise)
#       - total number of steps to perform  (optional, default value otherwise)



#call the initialization function from the initialization module



## While loop with multiple end criteria (nb of steps and total time  (possibly add convergence criteria) )
# 1. call the random number generator between 0 and 1(from numpy, random package, other?? in a function to change easily?) to get rho 
# 2. call the random number generator between 0 and 1(from numpy, random package, other??) to get xi
# 3. call the hop frequency calculator -> vector with frequency of possible events, sum of frequencies
# 4. call the call the event selection function -> selected event (number of neighbour which is swapped)
# 5. Accept event
#       a. update t = delta t
#       b. update step counter
#       c. update position array
#       d. update occupancy vector
#       e. update index array
# 6. Check if sampling at this step
#       a. if yes add current position array in next free layer of data collector

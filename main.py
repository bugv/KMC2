import global_functions
import json
import time

# Prompt user for input file
file_name = input("Please enter the name of the input file: ")

# load input structure from the json file
input_data = global_functions.read_input_file(file_name)
print("Input read from file")

# Run initialization
start_time = time.time()
init_struct = global_functions.initialization(*input_data)
end_time = time.time()
print("Total initialization Time", end_time - start_time)

# Output initialized structure
global_functions.write_full_to_json(init_struct, "Initialized_struct.json")
print("Initialized structure written to file")

# Run driver function

results = global_functions.driver(init_struct)
# Output results
start_time = time.time()
global_functions.write_full_to_json(results, "results.json")
end_time = time.time()
print("Driver function time", end_time - start_time)
# plot data???

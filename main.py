import global_functions
import json
import time
import psutil

# Ask user if create random alloy or used inputed structure
is_random_alloy = int(
    input(
        "1. The exact input structure should be used \n 2. A random alloy should be create from the composition"
    )
)

# Prompt user for input file
file_name = input("Please enter the name of the input file: ")


# load input structure from the json file (based on type of input)
print(is_random_alloy)
if is_random_alloy == 1:
    "The exact input structure will be used"
    input_data = global_functions.read_input_file(file_name)
if is_random_alloy == 2:
    input_data = global_functions.read_input_file_composition(file_name)
if is_random_alloy != 1 and is_random_alloy != 2:
    raise TypeError("Please select either 1 or 2 depending on the type of inputs")
print("Input read from file")

print(input_data)


# Run initialization
start_time = time.time()
process = psutil.Process()
mem_before = process.memory_info().rss / 1024  # in kilobytes
init_struct = global_functions.initialization(*input_data)
mem_after = process.memory_info().rss / 1024  # in kilobytes
mem_used = mem_after - mem_before
print("mem used for initialization (kilobytes)", mem_used)
end_time = time.time()
print("Total initialization Time", end_time - start_time)

# Output initialized structure
global_functions.write_full_to_json(init_struct, "Initialized_struct.json")
print("Initialized structure written to file")

# Run driver function

start_time = time.time()
process = psutil.Process()
mem_before = process.memory_info().rss / 1024  # in kilobytes
results = global_functions.driver(init_struct)
mem_after = process.memory_info().rss / 1024  # in kilobytes
mem_used = mem_after - mem_before
print("mem used for driver function(kilobytes)", mem_used)
end_time = time.time()
print("Driver function time", end_time - start_time)
# Output results
global_functions.write_full_to_json(results, "results.json")

# plot data???

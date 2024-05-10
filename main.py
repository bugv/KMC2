import global_functions
import json
import time
import psutil

# Prompt user for input file
file_name = input("Please enter the name of the input file: ")

# load input structure from the json file
input_data = global_functions.read_input_file(file_name)
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

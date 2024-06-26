import global_functions
import json
import time
import psutil
import argparse

total_time_start = time.time()


### Creating a better user interface input with a parser
parser = argparse.ArgumentParser(
    description="Process input json file and argument for type of file"
)
parser.add_argument("input_file", type=str, help="input file (json)")
parser.add_argument("-exact", action="store_true", help="Uses the exact input file")
parser.add_argument(
    "-binary", action="store_true", help="Creates a binary from the input file"
)

args = parser.parse_args()

if not args.binary and not args.exact:
    parser.error(
        "Use either the -exact or the -binary argument to indicate whether the input file contains an exact structure or a binary alloy to construct"
    )

if args.exact:
    input_data = global_functions.read_input_file(args.input_file)
if args.binary:
    input_data = global_functions.read_input_file_composition(args.input_file)

print(input_data)

### End parser user interface


### Initial user interface with promt

# # Ask user if create random alloy or used inputed structure
# is_random_alloy = int(
#     input(
#         "1. The exact input structure should be used \n 2. A random alloy should be create from the composition"
#     )
# )

# # Prompt user for input file
# file_name = input("Please enter the name of the input file: ")


# # load input structure from the json file (based on type of input)
# print(is_random_alloy)
# if is_random_alloy == 1:
#     "The exact input structure will be used"
#     input_data = global_functions.read_input_file(file_name)
# if is_random_alloy == 2:
#     input_data = global_functions.read_input_file_composition(file_name)
# if is_random_alloy != 1 and is_random_alloy != 2:
#     raise TypeError("Please select either 1 or 2 depending on the type of inputs")
# print("Input read from file")

# print(input_data)


### End Initial user interface with prompt

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

total_time_end = time.time()

print("The total time is of ", total_time_end - total_time_start)

# plot data???

""" 
This module executes the KMC.
It reads the input file, runs the initialization followed by the main loop and 
outputs a set of files containing the results initialized_struct.json (strcture 
and all arrays generated in the initialization), POSCAR file (to view the structure),
results.json (structure, and results after running the KMC)

There are no functions defined in this module.

"""

import global_functions
import time
import argparse
import datetime
import os

total_time_start = time.time()



### Creating a better user interface input with a parser
parser = argparse.ArgumentParser(
    description="Process input json file and argument for type of file"
)
parser.add_argument("input_file", type=str, help="input file (json)")
parser.add_argument("output", type=str, default = None,nargs ="?", help="Directory to save results. If None or not specified, results are only saved temporarily.")
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

# outputs
dt = datetime.datetime.now().isoformat(timespec="seconds").replace(":","-")
raw_results_directory = args.output

os.makedirs("tmp", exist_ok=True)
if args.output is not None:
    os.makedirs(raw_results_directory, exist_ok=True)


#print(input_data)


# Run initialization
start_time = time.time()

if args.output is None:
    poscar_path = "tmp/POSCAR"
else :
    poscar_path = os.path.join(raw_results_directory, dt + "-POSCAR-" + args.input_file.split("/")[-1].rstrip(".json"))

init_struct = global_functions.initialization(*input_data,poscar_path=poscar_path)
end_time = time.time()
print("Total initialization Time", end_time - start_time)

# Output initialized structure

# First save for analysis.py
global_functions.write_full_to_json(init_struct, "tmp/Initialized_struct.json")
# Second save for archiving
if args.output is not None:
    init_struct_filename = dt + "-initialized_struct-" + args.input_file.split("/")[-1]
    global_functions.write_full_to_json(init_struct, os.path.join(raw_results_directory, init_struct_filename))

print("Initialized structure written to file")

# Run driver function
start_time = time.time()
results = global_functions.driver(init_struct)
end_time = time.time()
print("Driver function time", end_time - start_time)

# Output results
# First save for analysis.py
global_functions.write_full_to_json(results, "tmp/results.json")
# Second save for archiving
if args.output is not None:
    results_filename = dt + "-results-" + args.input_file.split("/")[-1]
    global_functions.write_full_to_json(results, os.path.join(raw_results_directory, results_filename))

total_time_end = time.time()
print("The total time is of ", total_time_end - total_time_start)

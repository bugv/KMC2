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
import analysis_v2
import h5py
import numpy as np

total_time_start = time.time()



### Creating a better user interface input with a parser
parser = argparse.ArgumentParser(
    description="Process input json file and argument for type of file"
)
parser.add_argument("input_file", type=str, help="input file (json)")
parser.add_argument("output_analysis", type=str, default = None,nargs ="?", help="output file (h5)")
parser.add_argument("-output_intermediate", type=str, default = None,nargs ="?", help="Directory to save results. If None or not specified, intermediate results are not saved.")
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
os.makedirs("tmp", exist_ok=True)
if args.output_intermediate is not None:
    os.makedirs(raw_results_directory, exist_ok=True)

dt = datetime.datetime.now().isoformat(timespec="seconds").replace(":","-")
raw_results_directory = args.output_intermediate

output_analysis = args.output_analysis
if output_analysis is None:
    output_analysis = os.path.join("tmp", "results_first_processing.h5")


# Run initialization
start_time = time.time()

if raw_results_directory is not None:
    poscar_path = os.path.join(raw_results_directory, dt + "-POSCAR-" + args.input_file.split("/")[-1].rstrip(".json"))
else :
    poscar_path = None

init_struct = global_functions.initialization(*input_data,poscar_path=poscar_path)
end_time = time.time()
print("Total initialization Time", end_time - start_time)

# Save for archiving
if raw_results_directory is not None:
    init_struct_filename = dt + "-initialized_struct-" + args.input_file.split("/")[-1]
    global_functions.write_full_to_json(init_struct, os.path.join(raw_results_directory, init_struct_filename))

    print(f"Initialized structure written to file {init_struct_filename}")

# Run driver function
start_time = time.time()
results_KMC = global_functions.driver(init_struct)
end_time = time.time()
print("Driver function time", end_time - start_time)

# Save for archiving
if raw_results_directory is not None:
    results_filename = dt + "-results-" + args.input_file.split("/")[-1]
    global_functions.write_full_to_json(results_KMC, os.path.join(raw_results_directory, results_filename))

# analysis

tic = time.time()

data_collector = results_KMC["data_collector"]
time_collector = results_KMC["time_collector"]

occupancy_vector = results_KMC["atom_pos"].occupancy_vector
start_indices = results_KMC["atom_pos"].index_array

results_analysis = analysis_v2.single_calc(data_collector, occupancy_vector, start_indices, results_KMC["atom_key"])

toc = time.time()
print("Analysis time", toc - tic)

group_keys = ['frequency_vector', 'total_nb_steps', 'sampling_frequency', 'time','composition_dict', 'atom_key']

with h5py.File(output_analysis,"w") as h5file:
    for key, value in results_analysis.items():
        data_group = h5file.create_group(key)
        for subkey, subvalue in value.items():
            data_group.create_dataset(subkey, data=subvalue)

    h5file.create_dataset("time_collector", data=results_KMC["time_collector"][0]) #[0] because for some reason time_collector is a np.array([[data]]) situation

    group_input = h5file.create_group("input")
    for key in group_keys:
        value = results_KMC[key]
        if isinstance(value, dict):
            subgroup = group_input.create_group(key)
            for subkey, subvalue in value.items():
                subgroup.create_dataset(subkey, data=subvalue)
        else:
            group_input.create_dataset(key, data=value)


total_time_end = time.time()
print("The total time is of ", total_time_end - total_time_start)

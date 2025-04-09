""" 
Script to run the KMC on a set of structures with varying compositions. 
It creates the appropriate files, runs the KMC, calculates the L_ij parameters and plots the L_ij parameters as a function of the composition
"""

import json
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time



### Modify the values here

data = {
    "number_steps": 1000000,
    "sampling_frequency": 1000,
    "neighbour_radius": 2.9,
    "hop_frequencies": {"Al": 364990214.518247042, "Fe": 364990214.518247042},
    "lattice": [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]],
    "atoms": ["Fe", "Fe", "Fe", "Fe"],
    "coords": [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5]],
    "supercell_size": [5, 5, 5],
    "composition": {"Al": 0.1, "Fe": 0.9},
}

elem_i = "Al"
elem_j = "Fe"

lowest_comp = 0.0
highest_comp = 1.0
comp_step = 0.05

input_dir = "inputs"
tmp_raw_kmc_output_file = "tmp/results.json"
output_dir = "results_first_processing"

os.makedirs(output_dir, exist_ok=True)
os.makedirs(input_dir, exist_ok=True)

filename_list = []
output_dir_list = []

for al_composition in [0.5] : #np.arange(lowest_comp, highest_comp, comp_step) :
    fe_composition = 1.0 - al_composition
    data["composition"] = {"Al": al_composition, "Fe": fe_composition}

    input_filename = os.path.join(
        input_dir,
        f"composition_Al_{al_composition:.2f}_Fe_{fe_composition:.2f}.json",
    )

    output_dir_name = os.path.join(
        output_dir,
        f"composition_Al_{al_composition:.2f}_Fe_{fe_composition:.2f}",
    )
    os.makedirs(output_dir_name, exist_ok=True)

    with open(input_filename, "w") as f:
        json.dump(data, f, indent=4)

    print(f"Created {input_filename}")

    filename_list.append(input_filename)
    output_dir_list.append(output_dir_name)


for i,input_filename in enumerate(filename_list):
    output_dir = output_dir_list[i]
    for i in range(5000) :
        try:
            print("---")
            print(i) 

            command = ["python3", "main.py", "-binary", input_filename]

            subprocess.run(command)
            print(f"Running main.py on {input_filename}")
            
            tic = time.time()
            output_filename = os.path.join(output_dir, f"{i}.h5")
            subprocess.run(["python3", "analysis_v2.py", tmp_raw_kmc_output_file, output_filename])
            print(f"running analysis.py on {input_filename}")
            toc = time.time()
            print(f"Time taken for analysis: {toc - tic} seconds")

        except Exception as e:
            print(f"Error running  {input_filename}: {e}")


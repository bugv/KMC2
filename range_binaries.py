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

output_dir = "results"

### End modify values

os.makedirs(output_dir, exist_ok=True)

filename_list = []

al_composition = 0.0
for al_composition in np.arange(lowest_comp, highest_comp, comp_step):
    fe_composition = 1.0 - al_composition
    data["composition"] = {"Al": al_composition, "Fe": fe_composition}

    filename = os.path.join(
        "results",
        f"composition_Al_{al_composition:.2f}_Fe_{fe_composition:.2f}.json",
    )

    with open(filename, "w") as outfile:
        json.dump(data, outfile, indent=4)

    print(f"Created {filename}")

    filename_list.append(filename)
    al_composition += 0.05

print(filename_list)

for filename in filename_list:
    try:

        command = ["python3", "main.py", "-binary", filename]

        subprocess.run(command)
        print(f"Running main.py on {filename}")
        subprocess.run(["python3", "analysis.py"])
        print(f"running analysis.py on {filename}")
    except Exception as e:
        print(f"Error running  {filename}: {e}")


## Plot results
data = pd.read_csv("results/all_results.dat", delim_whitespace=True)

## format of this outputfile (created by analysis.py)
# name_elem0, name_elem1, name_elem2, comp_elem0, comp_elem1, comp_elem2, L_00, L_11, L01


check = data["name_elem0"] == elem_i

data.loc[
    check, ["comp_elem0", "comp_elem1", "name_elem0", "name_elem1", "L_00", "L_11"]
] = data.loc[
    check, ["comp_elem1", "comp_elem0", "name_elem1", "name_elem0", "L_11", "L_00"]
].values


plt.plot(data["comp_elem0"], data["L_11"], marker="o", label="L_AA")
plt.plot(data["comp_elem0"], data["L_00"], marker="o", label="L_BB")
plt.plot(data["comp_elem0"], data["L_01"], marker="o", label="L_AB")
plt.ylim(1e-12, 1e-6)
plt.xlabel("Concentration of B")
plt.ylabel("cm^2/s")
plt.yscale("log")
plt.legend()
plt.show()

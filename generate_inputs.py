import os
import json
import numpy as np

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

lowest_comp = 0.05
highest_comp = 1.0
comp_step = 0.05

input_dir = "inputs/GA1_GB10_binary"
tmp_raw_kmc_output_file = "tmp/results.json"
output_dir = "results_first_processing/GA1_GB10_binary"

os.makedirs(output_dir, exist_ok=True)
os.makedirs(input_dir, exist_ok=True)

filename_list = []
output_dir_list = []

for al_composition in np.arange(lowest_comp, highest_comp, comp_step) :
    fe_composition = 1.0 - al_composition
    data["composition"] = {"Al": al_composition, "Fe": fe_composition}
    data["hop_frequencies"]["Fe"] = data["hop_frequencies"]["Al"] * 10

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

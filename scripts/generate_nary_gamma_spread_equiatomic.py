import os
import json
import numpy as np

input_dir = f"inputs/gamma_spread"
hop_frequency = 364990214.518247042

data = {
    "number_steps": 1000000,
    "sampling_frequency": 1000,
    "neighbour_radius": 2.9,
    "hop_frequencies": {"Al": 364990214.518247042, "B": 364990214.518247042, "C": 364990214.518247042},
    "lattice": [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]],
    "atoms": ["Al", "Al", "Al", "Al"],
    "coords": [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5]],
    "supercell_size": [5, 5, 5],
    "composition": {"Al": 0.1, "Fe": 0.9},
}


os.makedirs(input_dir, exist_ok=True)



element_list = ["Al","B","C","Dy","Eu","F"]

for n in [3,4,5,6] :
    comp = 1/n
    for spread in np.linspace(1,3,5) :
        gammas = np.logspace(0,spread,n)
        data["composition"] = {e : comp for e in element_list[:n]}
        data["hop_frequencies"] = {e : hop_frequency * gammas[i] for i,e in enumerate(element_list[:n])}

        input_filename = os.path.join(
            input_dir,
            f"{n}ary_log_spread_{int(gammas[-1])}.json",
        )

        with open(input_filename, "w") as f:
            json.dump(data, f, indent=4)

        print(f"Created {input_filename}")


import os
import json
import numpy as np

input_dir = f"inputs/unary"
hop_frequency = 364990214.518247042

data = {
    "number_steps": 1000000,
    "sampling_frequency": 1000,
    "neighbour_radius": 2.9,
    "hop_frequencies": {"Al": 364990214.518247042},
    "lattice": [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]],
    "atoms": ["Al", "Al", "Al", "Al"],
    "coords": [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5]],
    "supercell_size": [5, 5, 5],
    "composition": {"Al": 1.0},
}


os.makedirs(input_dir, exist_ok=True)

gamma_list = []

for n in [3,4,5,6] :
    comp = 1/n
    for spread in np.linspace(1,3,5) :
        gammas = np.logspace(0,spread,n)

        for G in gammas :
            if G in gamma_list :
                continue

            gamma_list.append(G)
            #fcc 
            data["atoms"] = ["Al","Al","Al","Al"]
            data["coords"] = [[0,0,0],[0,0.5,0.5],[0.5,0.5,0],[0.5,0,0.5]]
            data["supercell_size"] = [5,5,5]
            data["neighbour_radius"] = 2.9
            data["hop_frequencies"] = {"Al":hop_frequency * G}

            input_filename = os.path.join(
                input_dir,
                f"unary_fcc_{G}.json",
            )

            with open(input_filename, "w") as f:
                json.dump(data, f, indent=4)

            #bcc 
            data["atoms"] = ["Al","Al"]
            data["coords"] = [[0,0,0],[0.5,0.5,0.5]]
            data["supercell_size"] = [6,6,7]
            data["neighbour_radius"] = 3.5
            data["hop_frequencies"] = {"Al":hop_frequency * G}
            input_filename = os.path.join(
                input_dir,
                f"unary_bcc_{G}.json",
            )
            with open(input_filename, "w") as f:
                json.dump(data, f, indent=4)

            #sc
            data["atoms"] = ["Al"]
            data["coords"] = [[0,0,0]]
            data["supercell_size"] = [8,8,8]
            data["neighbour_radius"] = 4.2
            data["hop_frequencies"] = {"Al":hop_frequency * G}
            input_filename = os.path.join(
                input_dir,
                f"unary_sc_{G}.json",
            )

        print(f"Created {input_filename}")


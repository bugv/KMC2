import os
import json
import numpy as np

input_dir = f"inputs/ideal_ternary"

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


X = (
    np.array([0.8, 0.7, 0.7, 0.6, 0.6, 0.6, 0.5, 0.5, 0.5, 0.5, 0.4, 0.4, 0.4,
        0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2,
        0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]),
    np.array([0.1, 0.2, 0.1, 0.3, 0.2, 0.1, 0.4, 0.3, 0.2, 0.1, 0.5, 0.4, 0.3,
            0.2, 0.1, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.7, 0.6, 0.5, 0.4, 0.3,
            0.2, 0.1, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]),
    np.array([0.1, 0.1, 0.2, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3,
            0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.1, 0.2, 0.3, 0.4, 0.5,
            0.6, 0.7, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
)

for xA,xB,xC in np.array(X).T :
    data["composition"] = {"Al": xA, "B": xB, "C": xC}

    input_filename = os.path.join(
        input_dir,
        f"composition_A_{xA:.2f}_B_{xB:.2f}_C_{xC:.2f}.json",
    )

    with open(input_filename, "w") as f:
        json.dump(data, f, indent=4)

    print(f"Created {input_filename}")


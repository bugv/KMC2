import os
import json
import numpy as np

original_dir = f"inputs/gamma_spread"
bcc_input_dir = f"inputs/gamma_spread_bcc"
sc_input_dir = f"inputs/gamma_spread_sc"

os.makedirs(bcc_input_dir, exist_ok=True)
os.makedirs(sc_input_dir, exist_ok=True)

for file in os.listdir(original_dir):
    if file.endswith(".json"):
        with open(os.path.join(original_dir, file), "r") as f:
            data = json.load(f)
        
        # bcc
        data["atoms"] = ["Al","Al"]
        data["coords"] = [[0,0,0],[0.5,0.5,0.5]]
        data["supercell_size"] = [6,6,7]
        data["neighbour_radius"] = 3.5

        new_filepath = os.path.join(bcc_input_dir, file)

        with open(new_filepath, "w") as f:
            json.dump(data, f, indent=4)
            
        # sc 
        data["atoms"] = ["Al"]
        data["coords"] = [[0,0,0]]
        data["supercell_size"] = [8,8,8]
        data["neighbour_radius"] = 4.2

        new_filepath = os.path.join(sc_input_dir, file)
        
        with open(new_filepath, "w") as f:
            json.dump(data, f, indent=4)


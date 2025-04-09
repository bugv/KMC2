import json
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import argparse
import global_functions

parser = argparse.ArgumentParser(
    description="Process KMC output json file into output h5 file"
)
parser.add_argument("input_file", type=str, help="input file (json)")
parser.add_argument("output", type=str, default = None,nargs ="?", help="output file (h5)")

input_file = parser.parse_args().input_file
output_file = parser.parse_args().output
if output_file is None:
    output_file = os.path.join("tmp", "results_first_processing.h5")

def single_calc(all_data_collector, indices_0, indices_1, indices_2):
    data_collector_0 = all_data_collector[:, indices_0, :]
    data_collector_1 = all_data_collector[:, indices_1, :]
    data_collector_2 = all_data_collector[:, indices_2, :]

    # Vectorized computation for vacancy
    displs_array_2 = data_collector_2[:, :, :] - data_collector_2[:, :, [0]]

    # Vectorized computation for Al & Fe
    displs_array_0 = data_collector_0[:, :, :] - data_collector_0[:, :, [0]]
    displs_array_1 = data_collector_1[:, :, :] - data_collector_1[:, :, [0]]

    # Reshape arrays in time-site-xyz order
    displs_array_2 = np.transpose(displs_array_2, (2, 1, 0))
    displs_array_0 = np.transpose(displs_array_0, (2,1,0))
    displs_array_1 = np.transpose(displs_array_1, (2,1,0))


    # sum over coordinates
    displs_per_time_0 = np.sum(displs_array_0, axis=1)
    displs_per_time_1 = np.sum(displs_array_1, axis=1)

    # dot product over coordinates for every time <sum_zeta,xi r_zeta(t) . r_xi(t)>
    theprod_00 = np.einsum('ij,ij->i', displs_per_time_0, displs_per_time_0)  # dot product over coordinates for every time
    theprod_01 = np.einsum('ij,ij->i', displs_per_time_0, displs_per_time_1)  # dot product over coordinates for every time
    theprod_11 = np.einsum('ij,ij->i', displs_per_time_1, displs_per_time_1) # dot product over coordinates for every time

    # r2s 
    r2_00 = np.einsum('ikj,ikj->ik', displs_array_0, displs_array_0)
    r2_11 = np.einsum('ikj,ikj->ik', displs_array_1, displs_array_1)
    r2_vacancy = np.einsum('ikj,ikj->ik', displs_array_2, displs_array_2)

    r2s_00 = np.sum(r2_00, axis=1)
    r2s_11 = np.sum(r2_11, axis=1)
    r2s_vacancy = np.sum(r2_vacancy, axis=1)

    # cross_terms
    cross_00 = theprod_00 - np.sum(r2_00, axis=1)
    cross_11 = theprod_11 - np.sum(r2_11, axis=1)


    return {"theprod_00": theprod_00, "theprod_01": theprod_01, "theprod_11": theprod_11, "r2s_00": r2s_00, "r2s_11": r2s_11, "r2s_vacancy": r2s_vacancy, "cross_00": cross_00, "cross_11": cross_11}

all_data = global_functions.read_full_from_json(input_file)
all_data_collector = all_data["data_collector"]
time_collector = all_data["time_collector"]

occupancy_vector = all_data["atom_pos"].occupancy_vector
start_indices = all_data["atom_pos"].index_array

indices_0 = start_indices[np.where(occupancy_vector == 0.0)[0]]
indices_1 = start_indices[np.where(occupancy_vector == 1.0)[0]]
indices_2 = start_indices[np.where(occupancy_vector == 2.0)[0]]

res = single_calc(all_data_collector, indices_0, indices_1, indices_2)

with open(input_file, "r") as i:
    input_data = json.load(i)

group_keys = ['frequency_vector', 'total_nb_steps', 'sampling_frequency', 'time','composition_dict', 'atom_key']
group = {key : all_data[key] for key in group_keys}

with h5py.File(output_file,"w") as h5file:
    for key, value in res.items():
        h5file.create_dataset(key, data=value)

    h5file.create_dataset("time_collector", data=all_data["time_collector"][0]) #[0] because for some reason time_collector is a np.array([[data]]) situation
    
    group_input = h5file.create_group("input")
    for key in group_keys:
        value = all_data[key]
        if isinstance(value, dict):
            subgroup = group_input.create_group(key)
            for subkey, subvalue in value.items():
                subgroup.create_dataset(subkey, data=subvalue)
        else:
            group_input.create_dataset(key, data=value)

    



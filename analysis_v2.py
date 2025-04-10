import json
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import argparse
import global_functions
import time

def single_calc(all_data_collector, occupancy_vector, start_indices,atom_key):

    displs_array = {}
    displs_per_time = {}
    r2s = {}
    for atom_str,int_label in atom_key.items(): 
        indices = start_indices[np.where(occupancy_vector == int_label)[0]]
        data_collector = all_data_collector[:, indices, :]
        displs_array[atom_str] = data_collector[:, :, :] - data_collector[:, :, [0]] 
        displs_array[atom_str] = np.transpose(displs_array[atom_str], (2, 1, 0)) # Reshape arrays in time-site-xyz order
        displs_per_time[atom_str] = np.sum(displs_array[atom_str], axis=1) # sum over sites

        r2 = np.einsum('ijk,ijk->ij', displs_array[atom_str], displs_array[atom_str])
        r2s[atom_str] = np.sum(r2, axis=1)

    theproduct = {}

    for atom_1 in displs_per_time.keys():
        for atom_2 in displs_per_time.keys():
            if atom_1 <= atom_2: #only one product per pair
                theproduct[atom_1 + atom_2] = np.einsum('ij,ij->i', displs_per_time[atom_1], displs_per_time[atom_2])


    return {"r2s": r2s, "theproduct": theproduct}

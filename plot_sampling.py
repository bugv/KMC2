import json
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import global_functions

def single_calc(all_data_collector, indices_0, indices_1, indices_2):
    data_collector_0 = all_data_collector[:, indices_0, :]
    data_collector_1 = all_data_collector[:, indices_1, :]
    data_collector_2 = all_data_collector[:, indices_2, :]

    #vacancy
    displs_array_2 = np.full((np.shape(data_collector_2)[2], np.shape(data_collector_2)[1],3), None)
    for time_step in range(data_collector_2.shape[2]):
        displs_array_2[time_step, 0, 0] = data_collector_2[0, 0, time_step] - data_collector_2[0, 0, 0]
        displs_array_2[time_step, 0, 1] = data_collector_2[1, 0,time_step] - data_collector_2[1, 0, 0]
        displs_array_2[time_step, 0, 2] = data_collector_2[2, 0, time_step] - data_collector_2[2, 0, 0]

    # Al & Fe
    displs_array_0 = np.full((np.shape(data_collector_0)[2], np.shape(data_collector_0)[1],3), None)
    displs_array_1 = np.full((np.shape(data_collector_1)[2], np.shape(data_collector_1)[1],3), None)

    for site in range(data_collector_0.shape[1]):
        for time_step in range(data_collector_0.shape[2]):  
            displs_array_0[time_step, site, 0] = data_collector_0[0, site, time_step] - data_collector_0[0, site, 0]
            displs_array_0[time_step, site, 1] = data_collector_0[1, site, time_step] - data_collector_0[1, site, 0]
            displs_array_0[time_step, site, 2] = data_collector_0[2, site, time_step] - data_collector_0[2, site, 0]

    for site in range(data_collector_1.shape[1]):
        for time_step in range(data_collector_1.shape[2]):
            displs_array_1[time_step, site, 0] = data_collector_1[0, site, time_step] - data_collector_1[0, site, 0]
            displs_array_1[time_step, site, 1] = data_collector_1[1, site, time_step] - data_collector_1[1, site, 0]
            displs_array_1[time_step, site, 2] = data_collector_1[2, site, time_step] - data_collector_1[2, site, 0]



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

    # cross_terms
    cross_00 = theprod_00 - np.sum(r2_00, axis=1)
    cross_11 = theprod_11 - np.sum(r2_11, axis=1)

    return {"theprod_00": theprod_00, "theprod_01": theprod_01, "theprod_11": theprod_11, "r2_00": r2_00, "r2_11": r2_11, "r2_vacancy": r2_vacancy, "cross_00": cross_00, "cross_11": cross_11}


r2s_00 = []
crosses_00 = []
r2s_11 = []
crosses_11 = []
theprod_01 = []
r2s_vacancy = []


for filename in os.listdir("raw_results"):
    if filename.endswith(".json") and "results" in filename and "Al_0.50_Fe_0.50_" in filename:
        print(filename)
        all_data = global_functions.read_full_from_json(os.path.join("raw_results", filename))
        all_data_collector = all_data["data_collector"]
        time_collector = all_data["time_collector"]

        occupancy_vector = all_data["atom_pos"].occupancy_vector
        start_indices = all_data["atom_pos"].index_array

        indices_0 = start_indices[np.where(occupancy_vector == 0.0)[0]]
        indices_1 = start_indices[np.where(occupancy_vector == 1.0)[0]]
        indices_2 = start_indices[np.where(occupancy_vector == 2.0)[0]]

        try : 
            res = single_calc(all_data_collector, indices_0, indices_1, indices_2)

            r2s_00.append(np.sum(res["r2_00"],axis=1))
            crosses_00.append(res["cross_00"])
            r2s_11.append(np.sum(res["r2_11"],axis=1))
            crosses_11.append(res["cross_11"])
            theprod_01.append(res["theprod_01"])
            r2s_vacancy.append(np.sum(res["r2_vacancy"],axis=1))

        except Exception as e:
            print(f"Error in {filename}: {e}")


r2s_00 = np.array(r2s_00,dtype=float)
r2s_11 = np.array(r2s_11,dtype=float)
crosses_00 = np.array(crosses_00,dtype=float)
total_00 = r2s_00 + crosses_00
theprod_01 = np.array(theprod_01,dtype=float)
crosses_11 = np.array(crosses_11,dtype=float)
total_11 = r2s_11 + crosses_11

#analytical expression
Gamma = 364990214.518247042 # Hz
a = 4 / 2**0.5 # AA
x_V = 0.002 
x_B = 0.5
x_A = 1-x_B
z = 12
d = 3 
rho = z / 2 / d
f = 0.7815

nat = 500

L_00 = x_V * x_A * rho * a**2 * Gamma * (1 - x_B / (1-x_V) *  (1-f))
predict_00 = L_00 * np.arange(1000) * 6 * 500 * 2.28023874e-04 / 1000

L_11 = x_V * x_B * rho * a**2 * Gamma * (1 - x_A / (1-x_V) *  (1-f))
L_01 = x_V * x_A * x_B / (1- x_V) * rho * a**2 * Gamma * (1-f)
predict_01 = L_01 * np.arange(1000) * 6 * 500 * 2.28023874e-04 / 1000
predict_11 = L_11 * np.arange(1000) * 6 * 500 * 2.28023874e-04 / 1000

# time 
plt.close()
plt.plot(time_collector[0])
plt.show()

# L00
plt.close()

plt.plot(r2s_00.mean(axis=0),label="mean R^2 terms (sum over all sites)")
plt.plot(crosses_00.mean(axis=0), label="mean cross terms")
plt.plot(total_00.mean(axis=0), label="mean total")

plt.fill_between(range(r2s_00.shape[1]), 
                 r2s_00.mean(axis=0) - r2s_00.std(axis=0), 
                 r2s_00.mean(axis=0) + r2s_00.std(axis=0), 
                 alpha=0.1)

plt.fill_between(range(crosses_00.shape[1]), 
                 crosses_00.mean(axis=0) - crosses_00.std(axis=0), 
                 crosses_00.mean(axis=0) + crosses_00.std(axis=0), 
                 alpha=0.3)

plt.plot(predict_00, label="predicted")
plt.plot(predict_01, label="predicted 01 cross term")
plt.plot(theprod_01.mean(axis=0), label="KMC 01 mean cross term")

plt.title(f"L_00 terms over {r2s_00.shape[0]} runs")

plt.legend()
plt.ylabel("AA^2")
plt.xlabel("propto time")
plt.grid()
plt.show()

# L11
plt.close()
plt.plot(r2s_11.mean(axis=0),label="mean R^2 terms (sum over all sites)")
plt.plot(crosses_11.mean(axis=0), label="mean cross terms")
plt.plot(total_11.mean(axis=0), label="mean total")

plt.fill_between(range(r2s_11.shape[1]), 
                 r2s_11.mean(axis=0) - r2s_11.std(axis=0), 
                 r2s_11.mean(axis=0) + r2s_11.std(axis=0), 
                 alpha=0.1)

plt.fill_between(range(crosses_11.shape[1]),
                    crosses_11.mean(axis=0) - crosses_11.std(axis=0),
                    crosses_11.mean(axis=0) + crosses_11.std(axis=0),
                    alpha=0.3)  


plt.plot(predict_11, label="predicted")
plt.title(f"L_11 terms over {r2s_11.shape[0]} runs")
plt.legend()
plt.grid()
plt.show()

#L10 

plt.close()
plt.plot(theprod_01.mean(axis=0),label="mean from KMC")
plt.fill_between(range(theprod_01.shape[1]), 
                 theprod_01.mean(axis=0) - theprod_01.std(axis=0), 
                 theprod_01.mean(axis=0) + theprod_01.std(axis=0), 
                 alpha=0.1)
plt.title(f"L01 terms over {theprod_01.shape[0]} runs")

plt.plot(predict_01, label="predicted")

plt.legend()
plt.grid()
plt.show()

# --- 
plt.close()
for i,r2_00 in enumerate(r2s_00[:10]) :
    print(r2_00)
    col = i%10
    print(col)
    plt.plot(r2_00,f"C{col}-")

for i,r2_11 in enumerate(r2s_11[:10]) :
    col = i%10
    plt.plot(r2_11,f"C{col}--")

plt.show()

## super long 

all_data = global_functions.read_full_from_json("raw_results/2025-03-21T00-36-26-results-Input_binary_superlong.json")
all_data_collector = all_data["data_collector"]
time_collector = all_data["time_collector"]

occupancy_vector = all_data["atom_pos"].occupancy_vector
start_indices = all_data["atom_pos"].index_array

indices_0 = start_indices[np.where(occupancy_vector == 0.0)[0]]
indices_1 = start_indices[np.where(occupancy_vector == 1.0)[0]]
indices_2 = start_indices[np.where(occupancy_vector == 2.0)[0]]

r2s_00_l = []
crosses_00_l = []
r2s_11_l = []
crosses_11_l = []
theprod_01_l = []
r2s_vacancy_l = []

for i in range(99) :
    print(i)
    time_slice = slice(i*1000,(i+1)*1000)
    print(time_slice)
    res = single_calc(all_data_collector[:,:,time_slice].copy(), indices_0, indices_1, indices_2)

    r2s_00_l.append(np.sum(res["r2_00"],axis=1))
    crosses_00_l.append(res["cross_00"])
    r2s_11_l.append(np.sum(res["r2_11"],axis=1))
    crosses_11_l.append(res["cross_11"])
    theprod_01_l.append(res["theprod_01"])
    r2s_vacancy_l.append(np.sum(res["r2_vacancy"],axis=1))

# Convert lists to numpy arrays
r2s_00_l = np.array(r2s_00_l, dtype=float)
r2s_11_l = np.array(r2s_11_l, dtype=float)
crosses_00_l = np.array(crosses_00_l, dtype=float)
total_00_l = r2s_00_l + crosses_00_l
theprod_01_l = np.array(theprod_01_l, dtype=float)
crosses_11_l = np.array(crosses_11_l, dtype=float)
total_11_l = r2s_11_l + crosses_11_l

# L00
plt.close()
plt.plot(r2s_00_l.mean(axis=0), label="mean R^2 terms (sum over all sites)")
plt.plot(crosses_00_l.mean(axis=0), label="mean cross terms")
plt.plot(total_00_l.mean(axis=0), label="mean total")

plt.fill_between(range(r2s_00_l.shape[1]),
                    r2s_00_l.mean(axis=0) - r2s_00_l.std(axis=0),
                    r2s_00_l.mean(axis=0) + r2s_00_l.std(axis=0),
                    alpha=0.1)

plt.fill_between(range(crosses_00_l.shape[1]),
                    crosses_00_l.mean(axis=0) - crosses_00_l.std(axis=0),
                    crosses_00_l.mean(axis=0) + crosses_00_l.std(axis=0),
                    alpha=0.3)

plt.plot(predict_00, label="predicted")
plt.title(f"L_00 terms over {r2s_00_l.shape[0]} runs")
plt.legend()
plt.ylabel("AA^2")
plt.xlabel("propto time")
plt.grid()
plt.show()

# L11
plt.close()
plt.plot(r2s_11_l.mean(axis=0), label="mean R^2 terms (sum over all sites)")
plt.plot(crosses_11_l.mean(axis=0), label="mean cross terms")
plt.plot(total_11_l.mean(axis=0), label="mean total")

plt.fill_between(range(r2s_11_l.shape[1]),
                    r2s_11_l.mean(axis=0) - r2s_11_l.std(axis=0),
                    r2s_11_l.mean(axis=0) + r2s_11_l.std(axis=0),
                    alpha=0.1)

plt.fill_between(range(crosses_11_l.shape[1]),
                    crosses_11_l.mean(axis=0) - crosses_11_l.std(axis=0),
                    crosses_11_l.mean(axis=0) + crosses_11_l.std(axis=0),
                    alpha=0.3)

plt.plot(predict_11, label="predicted")
plt.title(f"L_11 terms over {r2s_11_l.shape[0]} runs")
plt.legend()
plt.grid()
plt.show()

# L01
plt.close()
plt.plot(theprod_01_l.mean(axis=0), label="mean from KMC")
plt.fill_between(range(theprod_01_l.shape[1]),
                    theprod_01_l.mean(axis=0) - theprod_01_l.std(axis=0),
                    theprod_01_l.mean(axis=0) + theprod_01_l.std(axis=0),
                    alpha=0.1)

plt.plot(predict_01, label="predicted")
plt.title(f"L_01 terms over {theprod_01_l.shape[0]} runs")
plt.legend()
plt.grid()
plt.show()
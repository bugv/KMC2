"""Module to calculate the kinetic transport coefficients from the results
When run, the L_ij vales are output to the results.dat file

functions:
- calculate_L_ij
- L_ij_plot
"""

import global_functions
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def calculate_L_ij(
    atom_type_i: str,
    atom_type_j: str,
    data_i: np.array,
    data_j: np.array,
    time_collector: np.array,
    nb_sites: int,
    filename: str,
) -> float:
    nb_site_i = np.shape(data_i)[1]
    nb_site_j = np.shape(data_j)[1]
    displs_array_i = np.full((np.shape(data_i)[2], np.shape(data_i)[1]), None)
    displs_array_j = np.full((np.shape(data_j)[2], np.shape(data_j)[1]), None)
    data_j = data_j * 1e-8  # convert to cm
    data_i = data_i * 1e-8  # convert to cm

    for site in range(data_i.shape[1]):
        for time_step in range(data_i.shape[2]):
            displ = np.sqrt(
                (
                    (data_i[0, site, time_step] - data_i[0, site, 0]) ** 2
                    + (data_i[1, site, time_step] - data_i[1, site, 0]) ** 2
                    + (data_i[2, site, time_step] - data_i[2, site, 0]) ** 2
                )
            )
            displs_array_i[time_step, site] = displ

    for site in range(data_j.shape[1]):
        for time_step in range(data_j.shape[2]):
            displ = np.sqrt(
                (
                    (data_j[0, site, time_step] - data_j[0, site, 0]) ** 2
                    + (data_j[1, site, time_step] - data_j[1, site, 0]) ** 2
                    + (data_j[2, site, time_step] - data_j[2, site, 0]) ** 2
                )
            )
            displs_array_j[time_step, site] = displ

    displs_per_time_i = np.sum(displs_array_i, axis=1)
    displs_per_time_j = np.sum(displs_array_j, axis=1)

    mds = (displs_per_time_i * displs_per_time_j) / (nb_sites**2)
    mds = np.array(mds, dtype=float)
    linevalues = stats.linregress(time_collector[0, :], mds)
    slope, intercept, r_value, p_value, std_err = linevalues

    print(
        f"L_{atom_type_i}{atom_type_j}_tilde [cm2/s]",
        slope / (6),
        atom_type_i,
        atom_type_j,
    )

    L_ij_plot(time_collector, linevalues, mds, atom_type_i, atom_type_j, filename)
    return slope / 6


def L_ij_plot(
    timecollector: np.array,
    linevalues: tuple,
    mds: np.array,
    atomtype_i: str,
    atomtype_j: str,
    filename: str,
):
    slope, intercept, r_value, p_value, std_err = linevalues
    L_ij = slope / 6
    plt.clf()
    plt.plot(time_collector[0, :], mds[:], label="calculated")
    plt.plot(
        time_collector[0, :], intercept + slope * time_collector[0, :], label="fit"
    )
    plt.text(
        0.95,
        0.05,
        f"L_{atomtype_i}{atomtype_j} = {L_ij: .4e}",
        transform=plt.gca().transAxes,
        fontsize=12,
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox=dict(facecolor="white", alpha=0.5, edgecolor="black"),
    )
    plt.xlabel("time [s]")
    plt.ylabel("mean square displacement [cm^2]")
    plt.grid()
    plt.legend()
    # plt.savefig(filename)
    # plt.show()


all_data = global_functions.read_full_from_json("results.json")
all_data_collector = all_data["data_collector"]
time_collector = all_data["time_collector"]

occupancy_vector = all_data["atom_pos"].occupancy_vector

indices_0 = np.where(occupancy_vector == 0.0)[0]
indices_1 = np.where(occupancy_vector == 1.0)[0]
indices_2 = np.where(occupancy_vector == 2.0)[0]


nb_sites = np.shape(all_data_collector)[1]
print(np.shape(indices_0))
print(np.shape(indices_1))
print(np.shape(indices_2))


data_collector_0 = all_data_collector[:, indices_0, :]
data_collector_1 = all_data_collector[:, indices_1, :]
data_collector_2 = all_data_collector[:, indices_2, :]


for key, value in all_data["atom_key"].items():
    if value == 0:
        elem_0 = key
        comp_0 = all_data["composition_dict"][key]


for key, value in all_data["atom_key"].items():
    if value == 1:
        elem_1 = key
        comp_1 = all_data["composition_dict"][key]


for key, value in all_data["atom_key"].items():
    if value == 2:
        elem_2 = key
        comp_2 = all_data["composition_dict"][key]


L_00 = calculate_L_ij(
    elem_0,
    elem_0,
    data_collector_0,
    data_collector_0,
    time_collector,
    nb_sites,
    f"results/plot_L{elem_0}{elem_0}_comp_{elem_0}{comp_0}_{elem_1}{comp_1}.png",
)


L_11 = calculate_L_ij(
    elem_1,
    elem_1,
    data_collector_1,
    data_collector_1,
    time_collector,
    nb_sites,
    f"results/plot_L{elem_1}{elem_1}_comp_{elem_0}{comp_0}_{elem_1}{comp_1}.png",
)

L_01 = calculate_L_ij(
    elem_0,
    elem_1,
    data_collector_0,
    data_collector_1,
    time_collector,
    nb_sites,
    f"results/plot_L{elem_0}{elem_1}_comp_{elem_0}{comp_0}_{elem_1}{comp_1}.png",
)

with open("results/all_results.dat", "a") as file:
    file.write(
        f"{elem_0} {elem_1} {elem_2} {comp_0} {comp_1} {comp_2} {L_00} {L_11} {L_01}\n"
    )


# def get_Ls_as_paper(
#     atom_type_i: str,
#     atom_type_j: str,
#     data_i: np.array,
#     data_j: np.array,
#     time_collector: np.array,
#     nb_sites: int,
#     filename: str,
# ) -> float:
#     ### I have an array with for atoms of type i
#     # get the array of displacement vectors for each atom of type i at time t
#     # take the sum of the displacement of all the atoms at time t -> array with one displacement vector for each time step
#     # take the produit scalaire of this vector with the vector for atom type j at the same time step
#     # take the average of these values?
#     nb_site_i = np.shape(data_i)[1]
#     nb_site_j = np.shape(data_j)[1]
#     displs_array_i = np.full((np.shape(data_i)[2], np.shape(data_i)[1], 3), None)
#     displs_array_j = np.full((np.shape(data_j)[2], np.shape(data_j)[1], 3), None)
#     data_j = data_j * 1e-8  # convert to cm
#     data_i = data_i * 1e-8  # convert to cm

#     for site in range(data_i.shape[1]):
#         for time_step in range(data_i.shape[2]):
#             displ = (
#                 (data_i[0, site, time_step] - data_i[0, site, 0])
#                 + (data_i[1, site, time_step] - data_i[1, site, 0])
#                 + (data_i[2, site, time_step] - data_i[2, site, 0])
#             )
#             displs_array_i[time_step, site, :] = displ

#     for site in range(data_j.shape[1]):
#         for time_step in range(data_j.shape[2]):
#             displ = (
#                 (data_j[0, site, time_step] - data_j[0, site, 0])
#                 + (data_j[1, site, time_step] - data_j[1, site, 0])
#                 + (data_j[2, site, time_step] - data_j[2, site, 0])
#             )
#             displs_array_j[time_step, site, :] = displ

#     print(displs_array_i)
#     print(np.shape(displs_array_i))
#     displs_per_time_i = np.sum(displs_array_i, axis=0)
#     displs_per_time_j = np.sum(displs_array_j, axis=0)

#     # take the dot product for each time -> array with one value per time step

#     array_dot_prod = np.full((1, np.shape(data_i)[1]), None)
#     for t in range( np.shape(data_i)[1])
#         array_dot_prod[0, t] = displs_array_i[t,0]* displs_array_j[t,0] + displs_per_time_i[t,1]* displs_array_j[t,1] + displs_array_i[t,2]* displs_array_j[t,2]

#     print(array_dot_prod)

#     print(np.shape(displs_per_time_i))

# mds = (displs_per_time_i * displs_per_time_j) / (
#     nb_sites**2
# )  # should num sites be numsites_i +numsites_j
# mds = np.array(mds, dtype=float)
# # print("mds calculated", np.shape(mds))
# linevalues = stats.linregress(time_collector[0, :], mds)
# slope, intercept, r_value, p_value, std_err = linevalues
# # print("slope calculated", slope)
# # print(
# # "slope",
# # slope,
# # "intercept",
# # intercept,
# # "r value",
# # r_value,
# # "p value",
# # p_value,
# # "std err",
# # std_err,
# # )
# print(
#     f"L_{atom_type_i}{atom_type_j}_tilde [cm2/s]",
#     slope / (6),
#     atom_type_i,
#     atom_type_j,
# )

# L_ij_plot(time_collector, linevalues, mds, atom_type_i, atom_type_j, filename)
# return slope / 6

# # plt.plot(time_collector[0, :], mds[:], label="calculated")
# # plt.plot(
# #     time_collector[0, :], intercept + slope * time_collector[0, :], label="fit"
# # )
# # plt.xlabel("time [s]")
# # plt.ylabel("mean square displacement [cm^2]")
# # plt.grid()
# # plt.legend()
# # plt.show()


# get_Ls_as_paper(
#     elem_0,
#     elem_0,
#     data_collector_0,
#     data_collector_0,
#     time_collector,
#     nb_sites,
#     "testplot.png",
# )

# print("L_AlFe", L_ij)
##need to know which columns of the data collector contain the atoms of each type
## -> use the data occupancy vector


# def plot_mean_square_displ(filename: str) -> float:
#     data = global_functions.read_full_from_json(filename)
#     time_collector = data["time_collector"]
#     data_collector = data["data_collector"]
#     nb_sites = np.shape(data_collector)[1]
#     print("number sites", nb_sites)
#     displs_array = np.full(
#         (np.shape(data_collector)[2], np.shape(data_collector)[1]), None
#     )
#     nb_sites = np.shape(data_collector)[1]
#     data_collector = data_collector * 1e-8  # conver A to cm
#     for site in range(data_collector.shape[1]):
#         for time_step in range(data_collector.shape[2]):
#             displ = np.sqrt(
#                 (
#                     (data_collector[0, site, time_step] - data_collector[0, site, 0])
#                     ** 2
#                     + (data_collector[1, site, time_step] - data_collector[1, site, 0])
#                     ** 2
#                     + (data_collector[2, site, time_step] - data_collector[2, site, 0])
#                     ** 2
#                 )
#             )
#             displs_array[time_step, site] = displ
#     displs_array_without_vac = np.delete(
#         displs_array, data["initial_vac_position"], axis=1
#     )
#     mds = (np.sum(displs_array_without_vac, axis=1) ** 2) / (nb_sites**2)
#     mds = np.array(mds, dtype=float)
#     slope, intercept, r_value, p_value, std_err = stats.linregress(
#         time_collector[0, :], mds
#     )
#     print(
#         "slope",
#         slope,
#         "intercept",
#         intercept,
#         "r value",
#         r_value,
#         "p value",
#         p_value,
#         "std err",
#         std_err,
#     )
#     print(
#         "L_aa_tilde [cm2/s]",
#         slope / (6),
#     )

#     plt.plot(time_collector[0, :], mds[:], label="calculated")
#     plt.plot(
#         time_collector[0, :], intercept + slope * time_collector[0, :], label="fit"
#     )
#     plt.xlabel("time [s]")
#     plt.ylabel("mean square displacement (all atoms without vacancy) [cm^2]")
#     plt.grid()
#     plt.legend()
#     plt.show()


# plot_mean_square_displ("results.json")


#  # size = nb_time_steps x nb_sites
#     # each column corresponds to a site and each row a time step
#     # fill in this array, need to loop over each site and each time step
#     # need to repeat for each site at each time step -> loop over time steps
#     # then loop over the site in the structure
#     for site in range(data_collector.shape[1]):
#         for time_step in range(data_collector.shape[2]):
#             # get displ
#             displ_square = (
#                 (data_collector[0, site, time_step] - data_collector[0, site, 0]) ** 2
#                 + (data_collector[1, site, time_step] - data_collector[1, site, 0]) ** 2
#                 + (data_collector[2, site, time_step] - data_collector[2, site, 0]) ** 2
#             )
#             displs_array[time_step, site] = np.sqrt(
#                 displ_square
#             )  # have array with displacement of atom at each time step
#     displs_array_without_vac = np.delete(
#         displs_array, data["initial_vac_position"], axis=1
#     )
#     # for each row, take sum of displs
#     mds = np.sum(displs_array_without_vac, axis=1)
#     mds = (mds**2) / (nb_sites**2)
#     print("shapes, ", type(time_collector[0, :]), type(mds[:]))
#     print("shapes, ", np.shape(time_collector[0, :]), np.shape(mds[:]))
#     print("shapes, ", time_collector[0, :].dtype, mds[:].dtype)

#     new_time_collector = time_collector[0, :]
#     new_mds = np.array(mds, dtype=float)
#     slope, intercept, r_value, p_value, std_err = stats.linregress(
#         new_time_collector, new_mds
#     )
#     print(
#         "slope",
#         slope,
#         "intercept",
#         intercept,
#         "r value",
#         r_value,
#         "p value",
#         p_value,
#         "std err",
#         std_err,
#     )
#     print(
#         "L_aa_tilde [A2/s]",
#         slope / (6),
#     )

#     plt.plot(time_collector[0, :], mds[:], label="calculated")
#     plt.plot(
#         time_collector[0, :], intercept + slope * time_collector[0, :], label="fit"
#     )
#     plt.xlabel("time [s]")
#     plt.ylabel("mean square displacement (all atoms without vacancy) [Anstrom^2]")
#     plt.grid()
#     plt.legend()
#     plt.show()

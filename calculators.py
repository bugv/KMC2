## L_ij calculator
# s
import numpy as np
import math


##function to get the diffusion  constant Dv for a vacancy
# need t = time
# 1. get finial position of vacancy
# 2. get initial position of vacancy
# 3. get initial and final coordinates of the vacancy
# 4. get distance between initial and final coordinates of the vacancy
# 5. get final time
def get_diffusion_coeff(
    vac_pos: int, data_collector: np.array, time: float, nb_samples: int
) -> float:

    final_vac_coor = data_collector[:, vac_pos, nb_samples]
    initial_vac_coor = data_collector[:, vac_pos, 0]
    # distance = math.dist(final_vac_coor, initial_vac_coor)
    distance = initial_vac_coor - final_vac_coor
    distance = np.dot(distance, distance)
    d_v = distance**2 / (6 * time)
    print("final_vac_position", final_vac_coor)
    print("initial_vac_position", initial_vac_coor)
    print("distance", distance)
    print("time", time)
    return d_v

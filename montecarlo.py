import random

import numpy as np
from numba import njit

from grid import choose_random_direction, choose_random_monomer, illegal_move, move_monomer, get_cluster_grid
from grid_energy import calc_energy


@njit
def do_nothing(*args): pass


@njit
def thermal_fluctuation(new_energy, energy, beta):
    return random.uniform(0, 1) < np.exp(-beta * (new_energy - energy))


# 1 e)
@njit
def monte_carlo(grid, N_s, M, T, n=0, t_equil=np.inf, t_r=np.inf, on_iteration=do_nothing):
    epsilon = np.full(N_s, -69)
    energy = calc_energy(grid)
    epsilon[0] = energy
    beta = 1 / (T * 1.380649e-23)

    cluster_sizes = np.zeros(n)
    measure_index = 0

    for t in range(N_s):
        energy = calc_energy(grid)
        rand_monomer = choose_random_monomer(M)
        rand_direction = choose_random_direction()

        if illegal_move(grid, rand_monomer, rand_direction):
            continue
        else:
            new_grid = move_monomer(grid, rand_monomer, rand_direction)
            new_energy = calc_energy(new_grid)
            if new_energy < energy or thermal_fluctuation(new_energy, energy, beta):
                grid = new_grid
                energy = new_energy

        epsilon[t] = energy

        if (t - t_equil) > 0 and (t - t_equil) % t_r == 0:
            _, d = get_cluster_grid(grid)
            cluster_sizes[measure_index] = 2*M/d
            measure_index += 1

        if t % 50 == 0:
            on_iteration(grid, epsilon[:t])
    print("done!")
    return grid, epsilon, np.mean(cluster_sizes[:measure_index])

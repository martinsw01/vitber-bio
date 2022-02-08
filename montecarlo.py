import random

import numpy as np
from numba import njit

from grid import choose_random_direction, choose_random_monomer, illegal_move, move_monomer, get_cluster_grid
from grid_energy import calc_relative_energy, α


@njit
def do_nothing(*args): pass


@njit
def thermal_fluctuation(new_energy, energy, beta):
    return random.uniform(0, 1) < np.exp(-beta * (new_energy - energy))


@njit
def should_measure(t, t_equil, t_r):
    return (t - t_equil) > 0 and (t - t_equil) % t_r == 0


@njit
def choose_random_legal_move(grid, M):
    random_monomer = choose_random_monomer(M)
    random_direction = choose_random_direction()
    while illegal_move(grid, random_monomer, random_direction):
        random_monomer = choose_random_monomer(M)
        random_direction = choose_random_direction()
    return random_monomer, random_direction


# 1 e)
@njit
def monte_carlo(grid, N_s, M, T, n=0, t_equil=np.inf, t_r=np.inf, on_iteration=do_nothing):
    epsilon = np.full(N_s, -69)
    beta = 1 / (T * 1.380649e-23) * α / 2

    cluster_sizes = np.zeros(n)
    measure_index = 0

    for t in range(N_s):
        rel_energy = calc_relative_energy(grid)

        rand_monomer, rand_direction = choose_random_legal_move(grid, M)

        new_grid = move_monomer(grid, rand_monomer, rand_direction)
        new_rel_energy = calc_relative_energy(new_grid)

        if new_rel_energy < rel_energy or thermal_fluctuation(new_rel_energy, rel_energy, beta):
            grid = new_grid
            rel_energy = new_rel_energy

        epsilon[t] = rel_energy

        if should_measure(t, t_equil, t_r):
            _, d = get_cluster_grid(grid)
            cluster_sizes[measure_index] = 2 * M / d
            measure_index += 1

        if t % 100 == 1:
            on_iteration(grid, epsilon[:t])

    return grid, epsilon, np.mean(cluster_sizes[:measure_index])

import numpy as np
from grid import calc_energy, choose_random_polymer, choose_random_direction, is_illegal_move, move_polymer

k_b = 1


def thermal_fluctuations(new_energy, energy, T):
    beta = 1 / (k_b * T)
    return np.random.uniform(0, 1) < np.exp(-beta * (new_energy - energy))


def do_nothing(*args): pass


def monte_carlo(grid, grid_parameters, t_equil, t_r, n, on_iteration=do_nothing):
    """
    The Metropolis algorithm described here:
    https://www.math.ntnu.no/emner/TMA4320/2022v/prosjekter/Biofysikkprosjekt_2022.pdf
    Chapter 3.1 and 3.2, page 9 and 10

    :param grid_parameters: N, M, L - Grid size, polymers, and polymer size respectively
    :param t_equil: The number of time steps done before measurements of the systems starts.
    :param t_r: The number of time steps between each measurement.
    :param n: The number of measurements made after the equilibration of the system.
    :param on_iteration: Callback called each iteration (e.i. updating plot of grid and energy)
    :return: final grid, energy profile, mean grid after equilibrium
    """
    N_s = t_equil + t_r * n
    N, M, L = grid_parameters

    epsilon = np.empty(N_s)
    grids = np.empty((n, N, N))
    grids_index = 0

    for i in range(N_s):
        energy = calc_energy(grid)
        rand_polymer = choose_random_polymer(M)
        rand_direction = choose_random_direction()
        if is_illegal_move(grid, rand_polymer, rand_direction):
            continue
        else:
            new_grid = move_polymer(grid, rand_polymer, rand_direction)
            new_energy = calc_energy(grid)
            if new_energy < energy or thermal_fluctuations(new_energy, energy, T=273):
                energy = new_energy
                grid = new_grid
        if (i - t_equil) > 0 and (i - t_equil) % t_r == 0:
            grids[grids_index] = new_grid
            grids_index += 1
        epsilon[i] = energy
        on_iteration(grid, epsilon[:i])
    return grid, np.mean(grids, 0), epsilon


print("Martin er veldig, veldig kul egt")
print("Nesten like kul som Trygve og Philip!")

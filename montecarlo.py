import numpy as np
from grid_energy import calc_energy
import random
from grid import choose_random_direction, choose_random_monomer, illegal_move, move_monomer



#1 e)
def monte_carlo(grid, N_s, M, T):

    epsilon = np.zeros(N_s)
    energy = calc_energy(grid)
    epsilon[0] = energy
    beta = 1/(T*1.380649e-23)

    for i in range(N_s):
        energy = calc_energy(grid)
        rand_monomer = choose_random_monomer(M)
        rand_direction = choose_random_direction()

        if illegal_move(grid, rand_monomer, rand_direction):
            continue
        else:
            new_grid = move_monomer(grid, rand_monomer, rand_direction)
            new_energy = calc_energy(new_grid)
            if new_energy < energy or random.uniform(0, 1) < np.exp(-beta*(new_energy-energy)):
                grid = new_grid
                energy = new_energy

        epsilon[i] = energy

    return grid, epsilon




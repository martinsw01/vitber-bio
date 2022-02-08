import numpy as np

from grid import neighbor_coordinates
from numba import njit

ε_r = 78  # [1] Relative permittivity
ε_0 = 8.85418782e-12  # [m^-3 kg^-1 s^4 A^2] vacuum permittivity
a = 23e-6  # [m] Distance between grid points (23 microns)
e = 1.60217662e-19  # [C]
α = e ** 2 / (4 * np.pi * ε_r * ε_0 * a ** 2)


# 1 d)
@njit
def calc_energy(grid):
    N = len(grid)
    monomer_indices = np.where(grid != 0)
    relative_energy = 0

    for monomer_index in zip(*monomer_indices):
        monomer = grid[monomer_index]
        for neighbor_coord in neighbor_coordinates(N, *monomer_index):
            neighbour_monomer = grid[neighbor_coord]
            if neighbour_monomer != monomer and neighbour_monomer != 0:  # i.e. not the same type or solvent
                relative_energy += np.sign(monomer) * np.sign(neighbour_monomer)

    return relative_energy * α / 2

import numpy as np


def choose_random_polymer(M):
    return np.random.randint(-M, M)


directions = np.array([[-1, 0], [0, -1], [1, 0], [0, 1]])


def choose_random_direction():
    i = np.random.randint(0, 3)
    return directions[i]


def is_illegal_move(grid, polymer, direction):
    N = len(grid)
    polymer_indices = np.where(grid == polymer)
    new_polymer_indices = tuple(((np.array(polymer_indices).T + direction) % N).T)
    is_solvent = grid[new_polymer_indices] == 0
    is_same_polymer = grid[new_polymer_indices] == polymer
    return not np.logical_or(is_solvent, is_same_polymer).all()


def move_polymer(grid, polymer, direction):
    N = len(grid)
    new_grid = np.copy(grid)
    polymer_indices = np.where(new_grid == polymer)
    new_polymer_indices = tuple(((np.array(polymer_indices).T + direction) % N).T)
    new_grid[polymer_indices] = 0
    new_grid[new_polymer_indices] = polymer
    return new_grid


def calc_energy(grid):
    N = len(grid)
    monomer_indices = np.where(grid != 0)
    energy = 0

    for i, j in zip(*monomer_indices):
        monomer = grid[i, j]
        for dx, dy in directions:
            neighbour_monomer = grid[(i + dx) % N, (j + dy) % N]
            if neighbour_monomer != monomer and neighbour_monomer != 0:  # i.e. not the same type or solvent
                some_energy = np.sign(monomer) * np.sign(neighbour_monomer)
                energy += some_energy

    return energy

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
    return np.logical_or(is_solvent, is_same_polymer).all()


def move_polymer(grid, polymer, direction):
    N = len(grid)
    polymer_indices = np.where(grid == polymer)
    new_polymer_indices = tuple(((np.array(polymer_indices).T + direction) % N).T)
    grid[polymer_indices] = 0
    grid[new_polymer_indices] = polymer

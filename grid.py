import numpy as np


def neighbor_coordinates(N, i, j):
    neighbors = (((i - 1) % N, j), (i, (j - 1) % N), ((i + 1) % N, j), (i, (j - 1) % N))
    return neighbors


def choose_random_monomer(M):
    a = np.random.randint(-M, M)
    while not a:
        a = np.random.randint(-M, M)
    return a


directions = np.array([[-1, 0], [0, -1], [1, 0], [0, 1]])


def choose_random_direction():
    i = np.random.randint(0, 3)
    return directions[i]


def illegal_move(grid, monomer, direction):
    N = len(grid)
    monomer_location = np.where(grid == monomer)
    if grid[(monomer_location[0][0] + direction[0])%N, (monomer_location[1][0] + direction[1])%N]:
        return True
    return False


def move_monomer(grid, monomer, direction):
    N = len(grid)
    new_grid = np.copy(grid)
    monomer_location = np.where(grid == monomer)
    new_grid[monomer_location[0][0], monomer_location[1][0]] = 0
    new_grid[(monomer_location[0][0] + direction[0]) % N, (monomer_location[1][0] + direction[1]) % N] = monomer
    return new_grid


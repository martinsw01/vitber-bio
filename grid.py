import numpy as np
from numba import njit


@njit
def neighbor_coordinates(N, i, j):
    neighbors = (((i - 1) % N, j), (i, (j - 1) % N), ((i + 1) % N, j), (i, (j + 1) % N))
    return neighbors


@njit
def choose_random_monomer(M):
    a = np.random.randint(-M, M)
    while not a:
        a = np.random.randint(-M, M)
    return a


directions = np.array([[-1, 0], [0, -1], [1, 0], [0, 1]])


@njit
def choose_random_direction():
    i = np.random.randint(0, 3)
    return directions[i]


@njit
def illegal_move(grid, monomer, direction):
    N = len(grid)
    monomer_location = np.where(grid == monomer)
    return grid[(monomer_location[0][0] + direction[0]) % N, (monomer_location[1][0] + direction[1]) % N] != 0


@njit
def move_monomer(grid, monomer, direction):
    N = len(grid)
    new_grid = np.copy(grid)
    monomer_location = np.where(grid == monomer)
    new_grid[monomer_location[0][0], monomer_location[1][0]] = 0
    new_grid[(monomer_location[0][0] + direction[0]) % N, (monomer_location[1][0] + direction[1]) % N] = monomer
    return new_grid


# 1 g)
@njit
def get_cluster_grid(grid):
    cluster_grid = np.zeros_like(grid)
    monomer_positions = np.where(grid != 0)
    cluster = 1
    for monomer_position in zip(*monomer_positions):
        if cluster_grid[tuple(monomer_position)] == 0:
            cluster_grid[tuple(monomer_position)] = cluster
            add_neighbours_to_cluster(grid, cluster_grid, cluster, monomer_position)
            cluster += 1

    return cluster_grid, cluster


@njit
def add_neighbours_to_cluster(grid, cluster_grid, cluster, monomer_coordinate):
    N = len(grid)
    for neighbor_coordinate in neighbor_coordinates(N, *monomer_coordinate):
        if grid[neighbor_coordinate] != 0 and cluster_grid[neighbor_coordinate] == 0:
            cluster_grid[neighbor_coordinate] = cluster
            add_neighbours_to_cluster(grid, cluster_grid, cluster, neighbor_coordinate)

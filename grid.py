import numpy as np
from numba import njit


@njit
def neighbor_coordinates(N, i, j):
    neighbors = (((i - 1) % N, j), (i, (j - 1) % N), ((i + 1) % N, j), (i, (j + 1) % N))
    return neighbors


@njit
def choose_random_polymer(M):
    a = np.random.randint(-M, M+1)
    while not a:
        a = np.random.randint(-M, M+1)
    return a


directions = np.array([[-1, 0], [0, -1], [1, 0], [0, 1]])


@njit
def choose_random_direction():
    i = np.random.randint(0, 4)
    return directions[i]


# 1 g)
@njit
def get_cluster_grid(grid):
    cluster_grid = np.zeros_like(grid)
    monomer_positions = np.where(grid != 0)
    cluster = 0
    for monomer_position in zip(*monomer_positions):
        if cluster_grid[monomer_position] == 0:
            cluster += 1
            cluster_grid[monomer_position] = cluster
            add_neighbours_to_cluster(grid, cluster_grid, cluster, monomer_position)

    return cluster_grid, cluster


@njit
def add_neighbours_to_cluster(grid, cluster_grid, cluster, monomer_coordinate):
    N = len(grid)
    for neighbor_coordinate in neighbor_coordinates(N, *monomer_coordinate):
        if grid[neighbor_coordinate] != 0 and cluster_grid[neighbor_coordinate] == 0:
            cluster_grid[neighbor_coordinate] = cluster
            add_neighbours_to_cluster(grid, cluster_grid, cluster, neighbor_coordinate)


# 2 f)
@njit
def is_polymer_broken(grid, polymer):
    N = len(grid)
    monomer_positions = np.array([[i, j] for i in range(N) for j in range(N) if grid[i, j] == polymer])

    clean_grid = np.zeros((N, N))
    for monomer in zip(monomer_positions):
        clean_grid[monomer] = polymer
    new_grid = return_cluster(clean_grid, monomer_positions[0][0], monomer_positions[0][1], polymer)
    new_monomer_positions = np.array([[i, j] for i in range(N) for j in range(N) if new_grid[i, j] == -polymer])

    if len(monomer_positions) == len(new_monomer_positions):
        return False
    return True


"""for monomer in monomer_positions:
        is_monomer_alone = False
        neighbor_coord = neighbor_coordinates(N, monomer[0], monomer[1])
        for neighbor in neighbor_coord:
            if grid[neighbor]==polymer:
                is_monomer_alone = False
                break
            else:
                is_monomer_alone=True
        if is_monomer_alone:
            return True
    return False"""


@njit
def return_cluster(grid, i, j, polymer):
    N = len(grid)
    neighbours = neighbor_coordinates(N, i, j)

    for neighbour in neighbours:
        if grid[neighbour[0], neighbour[1]] == polymer:
            grid[neighbour[0], neighbour[1]] = -polymer
            grid = return_cluster(grid, neighbour[0], neighbour[1], polymer)
    return grid



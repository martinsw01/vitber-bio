import numpy as np
import random

from grid import neighbor_coordinates


def create_grid(N, M):
    grid = np.zeros((N, N), dtype=int)
    indices = [(i, j) for i in range(N) for j in range(N)]
    random_monomer_indices = random.sample(indices, 2 * M)
    monomers = [m for m in range(-M, M + 1) if m != 0]
    for monomer, index in zip(monomers, random_monomer_indices):
        grid[index] = monomer
    return grid


# 2 a)
def random_available_position_index(grid):
    indices = np.array(np.where(grid == 0)).T
    i = np.random.randint(len(indices))
    return tuple(indices[i])


def available_neighbor_positions(grid, monomer_positions):
    N = len(grid)
    return {position for monomer in monomer_positions
            for position in neighbor_coordinates(N, *monomer)
            if grid[position] == 0}


def random_available_neighbor_position(grid, monomer_positions):
    neighbor_coords = available_neighbor_positions(grid, monomer_positions)
    return random.sample(neighbor_coords, 1)[0]


def create_polymer_grid(N, M, L):
    grid = np.zeros((N, N), dtype=int)
    polymers = [polymers for polymers in range(-M, M + 1) if polymers != 0]

    for polymer in polymers:
        rand_position = random_available_position_index(grid)
        grid[rand_position] = polymer

        monomer_positions = np.zeros((L, 2), dtype=int)
        monomer_positions[0] = rand_position

        for j in range(1, L):
            try:
                neighbor_position = random_available_neighbor_position(grid, monomer_positions[:j])
                grid[neighbor_position] = polymer
                monomer_positions[j] = neighbor_position
            except ValueError as e:
                print(f"Failed to create an {N}x{N} grid with {M} polymers, each containing {L} monomers: ", e)
                print("Retrying ...")
                return create_polymer_grid(N, M, L)

    return grid

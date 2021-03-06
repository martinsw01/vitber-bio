import numpy as np
from numba import njit
from grid import is_polymer_broken


@njit
def is_legal_rigid_move(grid, polymer, attempted_new_positions):
    for i, j in attempted_new_positions:
        monomer = grid[i, j]
        if monomer != 0 and monomer != polymer:
            return False
    else:
        return True


@njit
def rigid_move(grid, polymer, direction):
    N = len(grid)
    monomer_positions = np.array([[i, j] for i in range(N) for j in range(N) if grid[i, j] == polymer])

    attempted_new_positions = (monomer_positions + direction) % N

    if is_legal_rigid_move(grid, polymer, attempted_new_positions):
        new_grid = np.copy(grid)
        for i, j in monomer_positions:
            new_grid[i, j] = 0
        for i, j in attempted_new_positions:
            new_grid[i, j] = polymer
        return new_grid
    else:
        return grid


# 2 e)
@njit
def medium_flexibility_move(grid, polymer, direction):
    N = len(grid)
    monomer_positions = np.array([[i, j] for i in range(N) for j in range(N) if grid[i, j] == polymer])
    attempted_new_positions = (monomer_positions + direction) % N
    if is_legal_rigid_move(grid, polymer, attempted_new_positions):
        return rigid_move(grid, polymer, direction)

    new_grid = np.copy(grid)
    for coordinate in monomer_positions:
        value_at_new_position = grid[(coordinate[0]+direction[0]) % N, (coordinate[1] + direction[1]) %N]

        if value_at_new_position != 0 and value_at_new_position != polymer:
            for ind, coord in enumerate(monomer_positions):
                if direction[0]:
                    if coord[1] == coordinate[1]:
                        attempted_new_positions[ind] -= direction
                else:
                    if coord[0] == coordinate[0]:
                        attempted_new_positions[ind] -= direction
    for i, j in monomer_positions:
        new_grid[i, j] = 0
    for i, j in attempted_new_positions:
        new_grid[i, j] = polymer
    check_grid = np.copy(new_grid)
    if is_polymer_broken(check_grid, polymer):
        return grid
    return new_grid


@njit
def move_monomer(grid, monomer, direction):
    N = len(grid)
    new_grid = np.copy(grid)
    monomer_location = np.where(grid == monomer)
    new_grid[monomer_location[0][0], monomer_location[1][0]] = 0
    new_grid[(monomer_location[0][0] + direction[0]) % N, (monomer_location[1][0] + direction[1]) % N] = monomer
    return new_grid


@njit
def illegal_move(grid, monomer, direction):
    N = len(grid)
    monomer_location = np.where(grid == monomer)
    return grid[(monomer_location[0][0] + direction[0]) % N, (monomer_location[1][0] + direction[1]) % N] != 0
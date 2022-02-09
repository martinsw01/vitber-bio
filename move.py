import numpy as np
from numba import njit


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


#2 e)
# @njit
def medium_flexibility_move(grid, polymer, direction):
    N = len(grid)
    monomer_positions = np.array([[i, j] for i in range(N) for j in range(N) if grid[i, j] == polymer])
    attempted_new_positions = (monomer_positions + direction) % N
    if is_legal_rigid_move(grid, polymer, attempted_new_positions):
        return rigid_move(grid, polymer, direction)

    moving_monomers = np.copy(monomer_positions)
    new_grid = np.copy(grid)
    for coordinate in monomer_positions:
        value_at_new_position = grid[(coordinate[0]+direction[0]) % N, (coordinate[1] + direction[1]) %N]
        if value_at_new_position != 0 and value_at_new_position != polymer:
            for coord in monomer_positions:
                if direction[0]:
                    if coord[1] == coordinate[1]:
                        # index_of_frozen_monomer = np.array([])
                        index_of_frozen_monomer = np.where(np.any(moving_monomers-coord, axis=1)==False)[0][0]
                        moving_monomers = np.delete(moving_monomers, index_of_frozen_monomer, axis=0)
                        # moving_monomers = remove_coordinate_from_array_of_coordinates(moving_monomers, coord)
                else:
                    if coord[0] == coordinate[0]:
                        index_of_frozen_monomer = np.where(np.any(moving_monomers - coord, axis=1) == False)[0][0]
                        moving_monomers = np.delete(moving_monomers, index_of_frozen_monomer, axis=0)
                        #moving_monomers = remove_coordinate_from_array_of_coordinates(moving_monomers, coord)
    for monomer in moving_monomers:
        new_grid[monomer[0], monomer[1]] = 0
    for monomer in moving_monomers:
        new_grid[(monomer[0] + direction[0]) % N, (monomer[1] + direction[1]) % N] = polymer
    return new_grid


"""@njit
def remove_coordinate_from_array_of_coordinates(array, coordinate):
    for ind, val in enumerate(array):
        if val == coordinate:
            return np.delete(array, ind, axis=0)"""
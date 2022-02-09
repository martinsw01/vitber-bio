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

def medium_flexibility_move(grid, polymer, direction):
    N = len(grid)
    monomer_positions = np.array(np.where(grid == polymer)).T
    attempted_rigid_move = tuple((monomer_positions + direction).T % N)
    if is_legal_rigid_move(grid, polymer, attempted_rigid_move):
        return rigid_move(grid, polymer, direction)
    print("Is not legal rigid move")
    frozen_monomers = np.array([[-1, -1]])
    new_grid = np.copy(grid)
    for coordinate in monomer_positions:
        print(coordinate)
        value_at_new_position = grid[(coordinate[0]+direction[0]) % N, (coordinate[1] + direction[1]) %N]
        if value_at_new_position != 0 and value_at_new_position != polymer:
            for coord in monomer_positions:
                if direction[0]:
                    if coord[1] == coordinate[1]:
                        frozen_monomers = np.append(frozen_monomers, [coord], axis=0)
                else:
                    if coord[0] == coordinate[0]:
                        frozen_monomers = np.append(frozen_monomers, [coord], axis=0)
    moving_monomers = np.array([])
    print(frozen_monomers)
    for coordinate in monomer_positions:
        print(coordinate)
        if coordinate not in frozen_monomers:
            print("Hei")
            new_grid[coordinate[0], coordinate[1]] = 0
            moving_monomers = np.append(moving_monomers, coordinate)
    for monomer in moving_monomers:
        new_grid[(monomer[0] + direction[0]) % N, (monomer[1] + direction[1]) % N] = polymer
    return new_grid


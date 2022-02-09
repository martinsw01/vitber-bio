import numpy as np


def is_legal_rigid_move(grid, polymer, attempted_new_positions):
    monomers_at_new_position = grid[attempted_new_positions]
    is_solvent = monomers_at_new_position == 0
    is_same_polymer = monomers_at_new_position == polymer
    return np.logical_or(is_solvent, is_same_polymer).all()


def rigid_move(grid, polymer, direction):
    monomer_positions = np.array(np.where(grid == polymer)).T
    attempted_new_positions = tuple((monomer_positions + direction).T % len(grid))
    if is_legal_rigid_move(grid, polymer, attempted_new_positions):
        new_grid = np.copy(grid)
        #print(new_grid[monomer_positions[:, 0], monomer_positions[:, 1]])
        """for i in monomer_positions:
            new_grid[i[0], i[1]] = 0"""
        new_grid[monomer_positions[:, 0], monomer_positions[:, 1]] = 0
        new_grid[attempted_new_positions] = polymer
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


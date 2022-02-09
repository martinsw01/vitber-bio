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
    frozen_monomers = np.array([])
    for coordinate in monomer_positions:
        print(coordinate+direction)
        if grid[coordinate+direction % N] != 0 and grid[coordinate+direction % N] != polymer:
            for coord in monomer_positions:
                if direction[0]:
                    if coord[0] == coordinate[0]:
                        frozen_monomers = np.append(frozen_monomers, coord)
                else:
                    if coord[1] == coordinate[1]:
                        frozen_monomers = np.append(frozen_monomers, coord)
                        


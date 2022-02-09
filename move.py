import numpy as np


def is_legal_rigid_move(grid, polymer, attempted_new_positions):
    monomers_at_new_position = grid[attempted_new_positions]
    is_solvent = monomers_at_new_position == 0
    is_same_polymer = monomers_at_new_position == polymer
    return np.logical_or(is_solvent, is_same_polymer).alltrue()


def rigid_move(grid, polymer, direction):
    monomer_positions = np.array(np.where(grid == polymer)).T
    attempted_new_positions = tuple((monomer_positions + direction).T % len(grid))
    if is_legal_rigid_move(grid, polymer, attempted_new_positions):
        new_grid = np.copy(grid)
        new_grid[monomer_positions] = 0
        new_grid[attempted_new_positions] = polymer
        return new_grid
    else:
        return grid

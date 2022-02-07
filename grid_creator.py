from random import random

import numpy as np


def create_grid(N, M):
    grid = np.zeros((N, N))
    indices = [(i, j) for i in range(N) for j in range(N)]
    random_monomer_indices = random.sample(indices, 2*M)
    monomers = [m for m in range(-M, M + 1) if m != 0]
    for monomer, index in zip(monomers, random_monomer_indices):
        grid[index] = monomer
    return grid

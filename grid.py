import numpy as np

def neighbor_coordinates(N, i, j):
    neighbors = (((i - 1) % N, j), (i, (j - 1) % N), ((i + 1) % N, j), (i, (j - 1) % N))
    return neighbors
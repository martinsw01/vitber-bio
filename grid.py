import numpy as np

def neighbor_coordinates(N):
    neighbors = np.empty([N, N, 4], tuple)
    for i in range(N):
        for j in range(N):
            neighbors[i, j, 0] = ((i - 1) % N, j)
            neighbors[i, j, 1] = (i, (j - 1) % N)
            neighbors[i, j, 2] = ((i + 1) % N, j)
            neighbors[i, j, 3] = (i, (j - 1) % N)
    return neighbors
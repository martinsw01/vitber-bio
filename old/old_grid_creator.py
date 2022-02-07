import matplotlib.pyplot as plt
import numpy as np

neighborhood_add_indices = np.array([[-1, 0], [0, -1], [1, 0], [0, 1]])


def generate_random_position(high):
    return np.random.randint(np.zeros_like(high), high, dtype=int)


def is_legal_polymer_placement(grid, position):
    return grid[position] == 0


def is_legal_monomer_placement(grid, position):
    return grid[position] == 0


def get_random_valid_position(grid, is_legal):
    N = len(grid)

    valid_positions = [(i, j) for i in range(N) for j in range(N) if is_legal(grid, (i, j))]

    i = np.random.randint(0, len(valid_positions))
    return valid_positions[i]


def random_element(lst, high=np.inf):
    i = np.random.randint(min(len(lst), high))
    return lst[i]


def get_random_neighbour_position(grid, monomer_position):
    N = len(grid)
    neighborhood = (monomer_position + neighborhood_add_indices) % N
    valid_positions = [(i, j) for i, j in neighborhood if grid[i, j] == 0]
    return random_element(valid_positions)


def random_available_position(grid):
    available_positions = np.where(grid == 0)
    i = np.random.randint(len(grid))
    return tuple(np.array(available_positions).T[i])


def random_available_neighbour(grid, monomer_positions):
    N = len(grid)
    neighbourhood = [((i + dx) % N, (j + dy) % N) for dx, dy in neighborhood_add_indices
                     for i, j in monomer_positions
                     if grid[(i + dx) % N, (j + dy) % N] == 0]
    return random_element(neighbourhood)


def create_grid(N, M, L):
    grid = np.zeros((N, N), dtype=int)
    polymers = np.arange(-M-1, M + 1)
    np.random.shuffle(polymers)
    for i in polymers:
        if i == 0:
            continue
        random_position = get_random_valid_position(grid, is_legal_polymer_placement)
        grid[random_position] = i
        monomer_positions = np.zeros((L, 2), dtype=int)
        monomer_positions[0] = random_position
        for j in range(1, L):
            try:
                random_neighbour = random_available_neighbour(grid, monomer_positions[:j])
                grid[random_neighbour] = i
                monomer_positions[j] = random_neighbour
            except ValueError as e:
                raise ValueError(f"Not enough space for 2*M=2*{M} polymers with L={L} monomers in a NxN={N}x{N} grid", e)
    return grid


if __name__ == '__main__':
    for _ in range(100):
        plt.pcolormesh(np.sign(create_grid(100, 10, 50)))
        plt.pause(1)
    plt.show()

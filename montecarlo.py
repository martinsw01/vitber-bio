import numpy as np

k_b = 1


def random(new_energy, energy, T):
    beta = 1/(k_b * T)
    return np.random.uniform(0, 1) < np.exp(-beta * (new_energy - energy))


def monte_carlo(grid, calc_energy, rand_polymer, rand_direction, is_illegal, move_polymer):
    """

    :param grid:
    :param calc_energy: Function calculating energy in grid. Returns a float
    :param rand_polymer: Function choosing random polymer from grid
    :param rand_direction: Function choosing a random direction in grid
    :param is_illegal: Function determining whether move is illegal
    :param move_polymer: Function moves polymer in grid
    :return ?:

    """
    N = len(grid)
    epsilon = np.zeros(N)
    energy = calc_energy(grid)
    epsilon[0] = energy

    for i in range(N):
        energy = calc_energy(grid)
        polymer = rand_polymer(grid)
        direction = rand_direction()
        if is_illegal(grid, polymer, direction):
            continue
        else:
            new_grid = move_polymer(grid, polymer, direction)
            new_energy = calc_energy(grid)
            if new_energy < energy or random(new_energy, energy, T=273):
                energy = new_energy
                grid = new_grid
        epsilon[i] = energy


print("Martin er veldig, veldig kul egt")
print("Nesten like kul som Trygve og Philip!")

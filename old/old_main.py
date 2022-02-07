import matplotlib.pyplot as plt

from old_montecarlo import monte_carlo
from old_grid_creator import create_grid
from visualization import create_plotter


def main():
    grid_parameters = 15, 25, 1
    grid = create_grid(*grid_parameters)
    # final_grid, mean_grid, energy_profile = monte_carlo(grid, grid_parameters, t_equil=50_000, t_r=10, n=10)
    final_grid, mean_grid, energy_profile = monte_carlo(grid, grid_parameters, t_equil=50_000, t_r=10, n=10,
                                                        on_iteration=create_plotter())
    _, (ax1, ax2, ax3) = plt.subplots(3)
    ax1.pcolormesh(final_grid)
    ax2.pcolormesh(mean_grid)
    ax3.plot(energy_profile)

    plt.show()


if __name__ == '__main__':
    main()

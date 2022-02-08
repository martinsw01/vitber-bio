from montecarlo import monte_carlo
from visualization import show_monomers, create_plotter
from grid_creator import create_grid
import time


# Oppgave 1f)
def simulate_monomers(N, M, N_s, T):
    grid = create_grid(N, M)
    grid, epsilon = monte_carlo(grid, N_s, M, T)
    show_monomers(grid, M)


# Remember to turn off plots in tool window in Preferences > Tools > Python Scientific
# to see the animation.
def simulate_monomers_with_visualizations(N, M, N_s, T):
    grid = create_grid(N, M)
    _, _ = monte_carlo(grid, N_s, M, T, create_plotter())


# Remember that old code is still stored in git. If you want to see the chanegs to this document,
# right click the this file's tab, then Git > Show History
def main():
    simulate_monomers_with_visualizations(15, 25, 50000, 500)


if __name__ == "__main__":
    main()

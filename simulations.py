import matplotlib.pyplot as plt
import numpy as np

from grid_creator import create_grid
from montecarlo import monte_carlo
from visualization import plot_energy_and_grid, create_plotter


# Oppgave 1f)
def simulate_monomers(N, M, N_s, T):
    grid = create_grid(N, M)
    grid, epsilon, _ = monte_carlo(grid, N_s, M, T)
    plot_energy_and_grid(epsilon, grid)


# Remember to turn off plots in tool window in Preferences > Tools > Python Scientific
# to see the animation.
def simulate_monomers_with_visualizations(N, M, N_s, T):
    grid = create_grid(N, M)
    _, _, _ = monte_carlo(grid, N_s, M, T, create_plotter())


def t_equil(T, t_max, s, T_l, C):
    return t_max * np.exp(-s * (T - T_l)) + C


# 1h)
def sim_mean_cluster_size():
    t_max = 100_000
    s = 1 / 200
    T_l = 100
    T_h = 1_000
    temperatures = np.linspace(T_l, T_h, 10)
    t_r = 1000
    C = 10_000
    t_equil(temperatures, t_max, s, T_l, C)
    N = 15
    M = 25
    n = 100
    mean_cluster_sizes = [monte_carlo(create_grid(N, M),
                                      int(t_equil(T, t_max, s, T_l, C) + t_r * n),
                                      M, T, n,
                                      int(t_equil(T, t_max, s, T_l, C)),
                                      t_r)[2]
                          for T in temperatures]
    plt.plot(temperatures, mean_cluster_sizes, linestyle="", marker="o")
    plt.show()

def main():
    T = [200, 500]
    N = 15
    M = 25
    N_s = 500
    # simulate_monomers_with_visualizations(15, 25, 5000, 500)
    simulate_monomers(N, M, N_s, T[0])


if __name__ == "__main__":
    sim_mean_cluster_size()

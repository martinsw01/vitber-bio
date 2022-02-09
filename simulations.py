import numpy as np
from numba import njit

from grid_creator import create_grid, create_polymer_grid
from montecarlo import monte_carlo
from move import rigid_move
from visualization import plot_energy_and_grid, create_plotter, plot_mean_cluster_sizes


# Oppgave 1f)
def simulate_monomers(N, M, N_s, T):
    grid = create_grid(N, M)
    grid, epsilon, _ = monte_carlo(grid, N_s, M, T)
    plot_energy_and_grid(epsilon, grid)


# Remember to turn off plots in tool window in Preferences > Tools > Python Scientific
# to see the animation.
def simulate_monomers_with_visualizations(N, M, N_s, T):
    grid = create_grid(N, M)
    _, _, _ = monte_carlo(grid, N_s, M, T, on_iteration=create_plotter())


def t_equil(T, t_max, s, T_l, C):
    return int(t_max * np.exp(-s * (T - T_l)) + C)


# 1h)
def sim_mean_cluster_size():
    t_max = 100_000
    s = 1 / 200
    T_l = 100
    T_h = 1_000
    temperatures = np.linspace(T_l, T_h, 10)
    t_r = 1000
    C = 10_000
    N = 15
    M = 25
    n = 100

    results = [monte_carlo(grid=create_grid(N, M),
                           N_s=t_equil(T, t_max, s, T_l, C) + t_r * n,
                           M=M, T=T, n=n,
                           t_equil=t_equil(T, t_max, s, T_l, C),
                           t_r=t_r)
               for T in temperatures]

    mean_cluster_sizes = [mean_cluster_size for _, _, mean_cluster_size in results]
    plot_mean_cluster_sizes(mean_cluster_sizes, temperatures)


@njit
def allways_false(*args): return False


# 2 d)
def simulation_with_polymers():
    N = 30
    L = 13
    M = 6
    T = 200
    N_s = 30_000
    grid = create_polymer_grid(N, M, L)
    final_grid, energy, _ = monte_carlo(grid, N_s, M, T, move_polymer=rigid_move, is_illegal_move=allways_false)
    plt.pcolormesh(final_grid[::-1, ])
    plt.show()


def main():
    T = [200, 500]
    N = 15
    M = 25
    N_s = 50_000
    # simulate_monomers_with_visualizations(M, N, N_s, T[0])
    simulate_monomers(N, M, N_s, T[0])


if __name__ == "__main__":
    sim_mean_cluster_size()

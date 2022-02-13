# import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
from numba import njit

from grid_creator import create_grid, create_polymer_grid
from montecarlo import monte_carlo
from move import rigid_move, medium_flexibility_move
from visualization import plot_energy_and_grid, create_plotter, plot_measurements


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

    cluster_sizes = [2 * M / measurements for _, _, measurements in results]
    mean_cluster_sizes = [np.mean(m) for m in cluster_sizes]
    std_devs = [np.std(m) for m in cluster_sizes]
    plot_measurements([(mean_cluster_sizes, std_devs)], temperatures)


@njit
def always_false(*args): return False


# 2 d)
def simulation_with_polymers():
    N = 30
    L = 13
    M = 6
    T = 200
    N_s = 30_000
    grid = create_polymer_grid(N, M, L)
    final_grid, energy, _ = monte_carlo(grid, N_s, M, T, move_polymer=rigid_move, is_illegal_move=always_false,
                                        on_iteration=create_plotter())
    plot_energy_and_grid(energy, final_grid)


# 2 g)
def simulation_with_polymers_using_medium_flexibility():
    N = 30
    L = 13
    M = 6
    T = 200
    N_s = 30_000
    grid = create_polymer_grid(N, M, L)
    final_grid, energy, _ = monte_carlo(grid, N_s, M, T,
                                        move_polymer=medium_flexibility_move,
                                        is_illegal_move=always_false)
    plot_energy_and_grid(energy, final_grid)


# 2 h)
def calculate_expected_values():
    T = 300
    t_r = 1_000
    N = 30
    M = 5
    n = 400
    t_eq = 5_000
    N_s = n * t_r + t_eq
    mean_cluster_sizes_per_L = np.zeros(13)
    mean_m = np.zeros(13)
    std_m = np.zeros(13)
    std_d_per_L = np.zeros(13)
    polymer_sizes = np.linspace(13, 39, 3, dtype=int)

    for i in range(1):
        print(i)
        grid = create_polymer_grid(N, M, polymer_sizes[i])
        _, energy, number_of_clusters = monte_carlo(grid, N_s, M, T, n=n, t_equil=t_eq, t_r=t_r,
                                                    move_polymer=medium_flexibility_move, is_illegal_move=always_false)
        np.savez(f"oppgave2h_L_{polymer_sizes[i]}.npz", energy=energy, number_of_clusters=number_of_clusters)
        mean_number_of_clusters = np.mean(number_of_clusters)
        monomers = 2 * M * polymer_sizes[i]
        cluster_sizes = monomers / number_of_clusters
        number_of_clusters_std = np.std(number_of_clusters)
        cluster_sizes_std = np.std(cluster_sizes)
        mean_cluster_sizes_per_L[i] = mean_number_of_clusters / polymer_sizes[i]
        mean_m[i] = mean_number_of_clusters
        std_m[i] = cluster_sizes_std
        std_d_per_L[i] = number_of_clusters_std
    plot_measurements([(mean_cluster_sizes_per_L, std_d_per_L), (mean_m, std_m)], polymer_sizes)


def load_2h_and_plot():
    mean_cluster_sizes_per_L_array = np.zeros(1)
    std_cluster_size_per_L_array = np.zeros(1)
    mean_number_of_clusters = np.zeros(1)
    std_number_of_clusters_array = np.zeros(1)
    polymer_sizes = np.linspace(13, 39, 3, dtype=int)[:1]
    M = 5

    for i in range(1):
        data = np.load(f"oppgave2h_L_{polymer_sizes[i]}.npz")
        energy = data["energy"]
        # plt.plot(energy)
        number_of_clusters = data["number_of_clusters"]
        plt.plot(number_of_clusters)
        mean_number_of_clusters = np.mean(number_of_clusters)
        monomers_in_grid = 2 * M * polymer_sizes[i]
        cluster_sizes = monomers_in_grid / number_of_clusters
        number_of_clusters_std = np.std(number_of_clusters)
        cluster_sizes_std = np.std(cluster_sizes)
        mean_cluster_sizes_per_L_array[i] = mean_number_of_clusters / polymer_sizes[i]
        mean_number_of_clusters[i] = mean_number_of_clusters
        std_cluster_size_per_L_array[i] = cluster_sizes_std
        std_number_of_clusters_array[i] = number_of_clusters_std
    plot_measurements([(mean_cluster_sizes_per_L_array, std_number_of_clusters_array),
                       (mean_number_of_clusters, std_cluster_size_per_L_array)], polymer_sizes)



def main():
    T = [200, 500]
    N = 15
    M = 25
    N_s = 500_000
    simulate_monomers_with_visualizations(M, N, N_s, T[0])
    # simulate_monomers(N, M, N_s, T[0])


if __name__ == "__main__":
    # simulation_with_polymers()
    # simulation_with_polymers_using_medium_flexibility()
    # calculate_expected_values()
    load_2h_and_plot()
    # check_std_deviation()
    # sim_mean_cluster_size()
    # main()
    # decide_number_of_samples()

    # n = np.array([5, 10, 15, 25, 50, 75, 100, 200, 300, 400, 500])
    # plt.figure(1)
    # for i in range(7, 12):
    #     data = np.load(f"find_n_{i}_samples.npz")
    #     plot()

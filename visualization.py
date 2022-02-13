import matplotlib.pyplot as plt


def create_axs():
    _, axs = plt.subplots(2)
    return axs


# Animation
def create_plotter(axs=None):
    if axs is None:
        axs = create_axs()
    ax_mesh, ax_plot = axs

    def update_plot(grid, energy):
        ax_mesh.clear()
        ax_plot.clear()
        # ax_mesh.pcolormesh(np.sign(grid[::-1,]))
        ax_mesh.pcolormesh(grid[::-1, ])
        ax_plot.plot(energy)
        plt.pause(0.001)

    return update_plot


# oppgave 1b)
def show_monomers(grid, M):
    """
    Input is the grid (grid) as well as the number of monomers (M).
    """
    plt.figure(1)
    plt.pcolormesh(grid[::-1, ])
    plt.title(f"{2 * M} monomers")
    plt.show()


def show_polymers(grid, num_of_polymers, title="Plot of polymers"):
    plt.figure(1)
    plt.title(title)
    plt.pcolormesh(grid[::-1, ] * 10)
    plt.title(f"{num_of_polymers} polymers")
    plt.show()


# Oppgave 1f)
def plot_energy_and_grid(energy, final_grid):
    ax1, ax2 = create_axs()
    ax1.pcolormesh(final_grid[::-1, ])
    ax2.plot(energy)
    plt.show()


def plot_mean_cluster_sizes(mean_cluster_sizes, std_devs, temperatures):
    plt.plot(temperatures, mean_cluster_sizes)
    plt.errorbar(temperatures, mean_cluster_sizes, std_devs)
    plt.show()


titles = ["Cluster sizes per polymer size", "Number of clusters"]


def plot_measurements(measurements, L):
    _, axs = plt.subplots(len(measurements))
    if not isinstance(axs, np.ndarray):
        axs = [axs]
    for (mean, std), title, ax in zip(measurements, titles, axs):
        ax.set_title(title)
        ax.set_xlabel("L")
        ax.plot(L, mean)
        ax.errorbar(L, mean, std, marker='^', color="C0")
    plt.tight_layout()
    plt.show()


def plot_two_energies(energy1, energy2, name_energy1, name_energy2):
    plt.figure(1)
    plt.title(f"{name_energy1} vs {name_energy2}")
    plt.plot(energy1, label=name_energy1)
    plt.plot(energy2, label=name_energy2)
    plt.legend()
    plt.show()
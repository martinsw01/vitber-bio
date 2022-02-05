import matplotlib.pyplot as plt
import numpy as np


def create_axs():
    _, axs = plt.subplots(2)
    return axs


def create_plotter(axs=None):
    if axs is None:
        axs = create_axs()
    ax_mesh, ax_plot = axs

    def update_plot(grid, energy):
        ax_mesh.clear()
        ax_plot.clear()
        ax_mesh.pcolormesh(np.sign(grid))
        ax_plot.plot(energy)
        plt.pause(0.001)

    return update_plot

# -*- coding: utf-8 -*-
import os
import numpy as np
import pyvista as pv
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.special import erf


def kelvin_to_fahrenheit(T):
    return (T - 273.15) * 9/5 + 32


def analytical(z, t, Ti=297.04, Ts=810.97, alpha=4.0e-07):
    T = Ts + (Ti - Ts) * erf(-z / np.sqrt(4 * alpha * t))
    return kelvin_to_fahrenheit(T)


def _get_profile(region):
    z = np.array(region.cell_centers().points[:, 2].tolist())
    T = np.array(region.cell_data["T"].tolist())
    return z, T


def get_profile(grid):
    block = grid["block"]["internalMesh"]
    media = grid["media"]["internalMesh"]

    z_block, T_block = _get_profile(block)
    z_media, T_media = _get_profile(media)

    z = np.concatenate((z_block, z_media))
    T = np.concatenate((T_block, T_media))
    return z, kelvin_to_fahrenheit(T)


def tabulate_profile(reader):
    matrix = []

    for this_time in reader.time_values[1:]:
        reader.set_active_time_value(this_time)
        matrix.append(get_profile(reader.read())[1])

    matrix = np.array(matrix)
    return matrix


def plot_evolution(z, times, data, z_max=0.05, saveas="evolution.png"):
    mask = z < z_max
    z_selected = z[mask]

    order = np.argsort(z_selected)
    z_selected = z_selected[order]

    selected = data[:, mask][:, order]
    times_s = np.array(times)
    times = times_s / 60

    z_lims = (z_selected[0], z_selected[-1])
    t_lims = (times[0], times[-1])
    extent = (z_lims[0], z_lims[1], t_lims[0], t_lims[1])

    Z, TT = np.meshgrid(z_selected, times_s)
    analytical_field = analytical(-Z, TT)

    plt.close("all")

    cb = plt.imshow(selected, aspect="auto", extent=extent,
                    origin="lower", cmap="hot")

    plt.contour(z_selected, times, selected, levels=[200.0],
                colors="cyan", linewidths=1.5)

    cs = plt.contour(z_selected, times, analytical_field, levels=[200.0],
                    colors="lime", linewidths=1.5, linestyles="dashed")

    plt.xlabel("z [m]")
    plt.ylabel("time [min]")
    plt.title("Temperature [°F]")

    plt.legend(handles=[
        Line2D([0], [0], color="cyan",  linewidth=1.5,
               label="Numerical 200 °F"),
        Line2D([0], [0], color="lime",  linewidth=1.5, linestyle="dashed",
               label="Analytical 200 °F"),
    ], loc=4, fontsize="small")

    cbar = plt.colorbar(cb, label="Temperature [°F]")
    cbar.set_ticks([100, 200, 400, 600, 800, 1000])

    plt.savefig(saveas, dpi=300)


def go_get_back(f):
    def wrapper(path):
        try:
            f(path)
        finally:
            os.chdir(os.path.dirname(__file__))
    return wrapper


@go_get_back
def post_process(wd):
    os.chdir(wd)
    reader = pv.POpenFOAMReader("case.foam")
    reader.set_active_time_value(reader.time_values[1])
    times = reader.time_values[1:]
    z, _ = get_profile(reader.read())
    data = tabulate_profile(reader)
    plot_evolution(z, times, data)


post_process("sandbox-model-1")
post_process("sandbox-model-2")
post_process("sandbox-model-3")

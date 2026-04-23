# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pyvista as pv
from majordome_utilities.plotting import plot2d
from matplotlib import pyplot as plt


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
    F = (T - 273.15) * 9/5 + 32
    return z, F


def tabulate_profile(reader):
    matrix = []

    for this_time in reader.time_values[1:]:
        reader.set_active_time_value(this_time)
        matrix.append(get_profile(reader.read())[1])

    matrix = np.array(matrix)
    return matrix


reader = pv.POpenFOAMReader("case.foam")

reader.set_active_time_value(reader.time_values[1])

z, _ = get_profile(reader.read())

data = tabulate_profile(reader)


selected = data[:, z < 0.05]

# plot = plot2d(z, data[-1, :])
# plot.axes[0].set_xlim(0, 0.2)
# plot.show()

z_lims = (0.00, 0.05)
t_lims = (reader.time_values[1] / 60, reader.time_values[-1] / 60)
extent = (z_lims[0], z_lims[1], t_lims[0], t_lims[1])

plt.close("all")
plt.imshow(selected, aspect="auto", extent=extent, origin="lower", cmap="hot")
plt.xlabel("z [m]")
plt.ylabel("time [min]")
plt.title("Temperature [°F]")
plt.colorbar(label="Temperature [°F]")
plt.show()
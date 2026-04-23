# -*- coding: utf-8 -*-
import pyvista as pv
import numpy as np
from majordome_utilities.plotting import plot2d


def get_profile(region):
    z = np.array(region.cell_centers().points[:, 2].tolist())
    T = np.array(region.cell_data["T"].tolist())
    return z, T


reader = pv.POpenFOAMReader("case.foam")
times = reader.time_values

this_time = times[-1]
reader.set_active_time_value(this_time)

grid = reader.read()

block = grid["block"]["internalMesh"]
media = grid["media"]["internalMesh"]

z_block, T_block = get_profile(block)
z_media, T_media = get_profile(media)

z = np.concatenate((z_block, z_media))
T = np.concatenate((T_block, T_media))
F = (T - 273.15) * 9/5 + 32

plot = plot2d(z, F)
plot.axes[0].set_xlim(0, 0.2)
# block.plot(scalars="T", cmap="hot", show_edges=True)
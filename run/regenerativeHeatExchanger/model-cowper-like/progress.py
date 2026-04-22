# -*- coding: utf-8 -*-
import os
import heat_recovery.cowper_like as cl
import matplotlib.pyplot as plt

plt.ion()
os.chdir("model-cowper-like")

# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------

origin = "0.000000e+00"
origin = "3.800000e-01"
origin = "2.000000e+00"

p = cl.plot_temperature("Initialization", origin, loc=3)

p = cl.plot_pressure("Initialization", origin)

p = cl.plot_flowrate("Initialization", origin, loc=3)

p = cl.plot_table("solid", "solidTemperature", "volFieldValue", time=origin)

# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------

import pyvista as pv
from heat_recovery.cowper_like import CowperLikePost
from heat_recovery.thermocline import ThermoclineModel

def plot_temperature_field(post):
    opts = post.get_options()
    pl = pv.Plotter()

    # pl.add_mesh(post.slice_fluid, scalars="T", cmap="hot", **opts)
    pl.add_mesh(post.slice_solid, scalars="T", cmap="hot", **opts)
    CowperLikePost.align_camera(pl, xc=0.025, zc=0.02, ps=0.017)

    pl.show()


model = ThermoclineModel("dimensioning.yaml")

post = CowperLikePost(scale=(1, 1, model.num_D_h / model.num_h_t))
post.load_state()

plot_temperature_field(post)
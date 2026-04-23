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


def plot_temperature_field(post, field):
    opts = post.get_options()
    pl = pv.Plotter()

    match field:
        case "fluid":
             pl.add_mesh(post.slice_fluid, scalars="T", cmap="hot", **opts)
        case "solid":
             pl.add_mesh(post.slice_solid, scalars="T", cmap="hot", **opts)
        case _:
            pl.add_mesh(post.slice_fluid, scalars="T", cmap="hot", **opts)
            pl.add_mesh(post.slice_solid, scalars="T", cmap="hot", **opts)

    CowperLikePost.align_camera(pl, xc=0.025, zc=0.02, ps=0.017)
    pl.show()


def plot_field_by_key(post, name):
    opts = post.get_options()
    pl = pv.Plotter()
    pl.add_mesh(post.slice_fluid, scalars=name, cmap="jet", **opts)
    CowperLikePost.align_camera(pl, xc=0.0125, zc=0.02, ps=0.017)
    pl.show()


model = ThermoclineModel("dimensioning.yaml")

stretch_z = model.num_D_h / model.num_h_t

post = CowperLikePost(scale=(1, 1, stretch_z))
post.load_state()

plot_temperature_field(post, "solid_fluid")

plot_field_by_key(post, "p")
plot_field_by_key(post, "p_rgh")
plot_field_by_key(post, "rho")

fluid_base = post._fluid_internal.copy()
U_original = fluid_base["U"].copy()

fluid = fluid_base.scale((1, 1, stretch_z), inplace=False)
fluid = fluid.sample(fluid_base)

# fluid.cell_data["U_original"] = U_original

solid = post._solid_internal.copy()
solid = solid.scale((1, 1, stretch_z), inplace=False)

fluid.plot(scalars="U", cmap="jet", component=2)

# fluid.plot(scalars="T")
# solid.plot(scalars="T", cmap="hot")

# opts = post.get_options()
# pl = pv.Plotter()
# pl.add_mesh(post.slice_fluid, scalars="Umag", cmap="jet", **opts)
# pl.add_mesh(post.slice_fluid, scalars="U", component=2, cmap="jet", **opts)
# pl.add_mesh(post.slice_fluid, scalars="U", component=2, cmap="jet", **opts)
# CowperLikePost.align_camera(pl, xc=0.0125, zc=0.02, ps=0.017)
# pl.show()

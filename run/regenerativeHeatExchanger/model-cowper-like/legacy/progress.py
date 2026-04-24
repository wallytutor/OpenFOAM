# -*- coding: utf-8 -*-
import os
import heat_recovery.cowper_like as cl
import matplotlib.pyplot as plt

plt.ion()
os.chdir("model-cowper-like")

# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------

origin = "0.00000000e+00"
# origin = "3.800000e-01"
# origin = "2.000000e+00"

p = cl.plot_temperature("Initialization", origin, loc=3)

p = cl.plot_pressure("Initialization", origin)

p = cl.plot_flowrate("Initialization", origin, loc=3)

p = cl.plot_table("solid", "solidTemperature", "volFieldValue", time=origin)

p = cl.plot_table("solid", "solidEnergy", "volFieldValue", time=origin)

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



from pathlib import Path


class DomainPostprocessing:
    def __init__(self, domain):
        self.domain = domain
        self.reports = self._get_domain_reports()

    @property
    def domain_dir(self):
        if not hasattr(self, "_domain_dir"):
            self._domain_dir = Path("postProcessing") / self.domain
        return self._domain_dir

    def _get_domain_reports(self):
        if not self.domain_dir.is_dir():
            raise ValueError(
                f"Domain '{self.domain}' not found in postProcessing directory."
            )

        return [d.name for d in self.domain_dir.iterdir() if d.is_dir()]

    # def load_report(self, report):
    #     report_dir = self.domain_dir / report

    #     if not report_dir.is_dir():
    #         raise ValueError(f"Report '{report}' not found for domain '{self.domain}'.")

    #     # Implement loading logic here (e.g., read data files, parse results, etc.)
    #     print(f"Loading report '{report}' for domain '{self.domain}' from {report_dir}")
    #     # Placeholder: return the path to the report directory
    #     return report_dir




post = DomainPostprocessing("fluid")

# report = "flowRateInlet"
# post = Path("postProcessing")

# root_dir = post / domain / report
# sub_dirs = [d for d in root_dir.iterdir() if d.is_dir()]
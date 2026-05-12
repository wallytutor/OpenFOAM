# -*- coding: utf-8 -*-
""" Case automated post-processing.

**Important**: the boundary *outlet* is actually used as inlet and
vice-versa. This is because the case is set up for the charging phase,
where the flow enters through the outlet and exits through the *inlet*.

- [x] Inlet mass flow rate
- [x] Mass imbalance
- [x] Pressure drop
- [ ] Pressure profile
- [ ] Velocity profile
- [ ] Temperature profile

"""
import numpy as np
import pyvista as pv
from heat_recovery.thermocline import ThermoclineModel
from majordome import constants
from majordome.simulation import FoamPostProcessingLoader
from majordome.utilities import plot_xy
from pathlib import Path

RESULTS_DIR = Path("postProcessing")


def plot_inlet_mass_flow_rate(post):
    """ Plot the inlet mass flow rate over time. """
    df1 = post.load_report("flowRateOutlet")
    df2 = post.load_report("flowRateInlet")

    n = min(len(df1), len(df2))

    x = df1.iloc[:n, 0]
    m1 = df1.iloc[:n, 1] * (-1000)
    m2 = df2.iloc[:n, 1] * (-1000)

    p = plot_xy()
    p.axes[0].plot(x, m1, label="Outlet")
    p.axes[0].plot(x, m2, label="Inlet")
    p.axes[0].legend(loc="best")
    p.axes[0].set_xlabel("Time [s]")
    p.axes[0].set_ylabel("Mass flow rate [g/s]")
    p.savefig(RESULTS_DIR / "01-inlet_mass_flow_rate.png")


def plot_imbalance_mass_flow_rate(post):
    """ Plot the imbalance of mass flow rate over time. """
    df1 = post.load_report("flowRateOutlet")
    df2 = post.load_report("flowRateInlet")

    n = min(len(df1), len(df2))

    x = df1.iloc[:n, 0]
    m1 = df1.iloc[:n, 1] * (-1000)
    m2 = df2.iloc[:n, 1] * (-1000)
    imbalance = (m1 + m2) / m1 * 100

    p = plot_xy(x, imbalance)
    p.axes[0].set_xlabel("Time [s]")
    p.axes[0].set_ylabel("Mass flow imbalance [%]")
    p.savefig(RESULTS_DIR / "02-imbalance_mass_flow_rate.png")


def plot_total_pressure(post):
    """ Plot the total pressure over time. """
    df1 = post.load_report("pressureInlet")
    df2 = post.load_report("pressureOutlet")

    n = min(len(df1), len(df2))

    x = df1.iloc[:n, 0]
    p1 = (df1.iloc[:n, 1] - constants.P_NORMAL) / 100
    p2 = (df2.iloc[:n, 1] - constants.P_NORMAL) / 100

    p = plot_xy()
    p.axes[0].plot(x, p1, label="Outlet")
    p.axes[0].plot(x, p2, label="Inlet")
    p.axes[0].legend(loc="best")
    p.axes[0].set_xlabel("Time [s]")
    p.axes[0].set_ylabel("Total pressure [mbar]")
    p.savefig(RESULTS_DIR / "03-total_pressure.png")


def plot_pressure_drop(post):
    """ Plot the pressure drop over time. """
    df1 = post.load_report("pressureInlet")
    df2 = post.load_report("pressureOutlet")

    n = min(len(df1), len(df2))

    x = df1.iloc[:n, 0]
    p1 = (df1.iloc[:n, 1] - constants.P_NORMAL) / 100
    p2 = (df2.iloc[:n, 1] - constants.P_NORMAL) / 100
    dp = p2 - p1

    p = plot_xy(x, dp)
    p.axes[0].set_xlabel("Time [s]")
    p.axes[0].set_ylabel("Pressure drop [mbar]")
    p.savefig(RESULTS_DIR / "04-pressure_drop.png")


def get_options(time, scalar_bar={}, **kws) -> dict:
    fn_title = kws.pop("fn_title", lambda t: f"At t = {t:.0f} s\n")

    scalar_bar_args = {
        "title": fn_title(time),
        "vertical": False,
        "title_font_size": 16,
        "label_font_size": 12,
        "position_x": 0.1,
        "position_y": 0.05,
        "color": "k",
        "n_colors": 10,
        "width": 0.8,
        "height": 0.1,
        "fmt": "%.0f"
    }
    scalar_bar_args.update(scalar_bar)

    opts = {
        "pbr": False,
        "show_edges": False,
        "scalar_bar_args": scalar_bar_args,
    }
    opts.update(kws)

    return opts


def camera_position_xz(w: float, h: float) -> tuple:
    """ Get camera position for XZ plane view (top-down). """
    x = w / 2
    z = h / 2

    # Camera position above the geometry
    pos = (x, 2 * max(w, h), z)

    # Looking at the center
    focal = (x, 0, z)

    # View up: (0, 1, 0) - Y axis points up
    viewup = (0, 0, 1)

    return (pos, focal, viewup)


def align_camera(plot, xc=0.025, zc=0.02, ps=0.017):
    plot.renderer.set_background("#FFFFFF")
    plot.camera_position = camera_position_xz(xc, zc)
    plot.camera.parallel_projection = True
    plot.camera.parallel_scale = ps


def load_prepare_fields():
    model = ThermoclineModel("../dimensioning.yaml")
    reader = pv.POpenFOAMReader("case.foam")
    scale = (1, 1, model.num_D_h / model.num_h_t)

    time = reader.time_values[-1]
    reader.set_active_time_value(time)

    mesh = reader.read()

    data = mesh.slice("y")[0]
    data.points *= scale

    data.cell_data["rhoG"] = data.cell_data["rho"] * 1000.0
    data.cell_data["pRel"] = data.cell_data["p"] - constants.P_NORMAL
    data = data.cell_data_to_point_data()
    return data, time


def plot_field_pressure(data, time):
    plot = pv.Plotter()
    opts = get_options(time, scalar_bar={"fmt": "%.0f"})
    plot.add_mesh(data, scalars="pRel", cmap="jet", **opts)
    align_camera(plot, xc=0.0125, zc=0.02, ps=0.017)
    plot.show()


def plot_field_density(data, time):
    plot = pv.Plotter()
    opts = get_options(time, scalar_bar={"fmt": "%.2f"})
    plot.add_mesh(data, scalars="rhoG", cmap="jet", **opts)
    align_camera(plot, xc=0.0125, zc=0.02, ps=0.017)
    plot.show()


def plot_field_velocity(data, time):
    plot = pv.Plotter()
    opts = get_options(time, scalar_bar={"fmt": "%.2f"})
    plot.add_mesh(data, scalars="U", cmap="jet", **opts)
    align_camera(plot, xc=0.0125, zc=0.02, ps=0.017)
    plot.show()


if __name__ == "__main__":
    post = FoamPostProcessingLoader()
    plot_inlet_mass_flow_rate(post)
    plot_imbalance_mass_flow_rate(post)
    plot_total_pressure(post)
    plot_pressure_drop(post)

    fields, time = load_prepare_fields()
    plot_field_pressure(fields, time)
    plot_field_density(fields, time)
    plot_field_velocity(fields, time)

# -*- coding: utf-8 -*-

import os

# Mute C-level stderr (File Descriptor 2) immediately
_devnull = os.open(os.devnull, os.O_WRONLY)
os.dup2(_devnull, 2)
os.close(_devnull)

# Force headless off-screen environment
os.environ["VTK_DEFAULT_RENDER_WINDOW_HEADLESS"] = "1"
os.environ["PYVISTA_OFF_SCREEN"] = "true"
os.environ["PYVISTA_USE_IPYVISTA"] = "false"

import cantera as ct
import majordome as mj
import numpy as np
import pyvista as pv

from pathlib import Path

from numpy.typing import NDArray
from ruamel.yaml import YAML
from tabulate import tabulate

pv.set_jupyter_backend("static")

# TODO automate reloading if required:
_MESH = None

_SETUP = mj.FoamDictFile(path="model/0/_shared")


def load_yaml(fname: str | Path) -> dict:
    """ Load a YAML file. """
    with open(fname, encoding="utf-8") as f:
        config = YAML().load(f)
    return config


def get_power_supply(
        mech: ct.Solution,
        qdot_fuel: float,
        X_fuel: str | dict[str, float],
        X_oxid: str | dict[str, float],
        phase: str,
    ) -> mj.CombustionPowerSupply:
    """ Evaluate equivalent combustion power supply. """
    normal_flow = mj.NormalFlowRate(mech, X=X_fuel, name=phase)
    mdot_fuel = normal_flow(qdot_fuel) / 4

    chon = mj.CombustionAtmosphereCHON(mech, phase=phase, basis="mole")
    lhv = chon.solution_heating_value(X_fuel, X_oxid)

    supply = mj.CombustionPowerSupply(
        power       = lhv * mdot_fuel * 1000,
        equivalence = 1.0,
        fuel        = X_fuel,
        oxidizer    = X_oxid,
        mechanism   = mech,
        phase       = phase,
        basis       = "mole",
    )
    return supply


def load_setup_shared(
        air_temp: float = 0.0,
        fname: str | Path = "initialize.yaml",
        save: bool = False,
    ):
    """ Load simulation setup for the shared files. """
    conf = load_yaml(fname)

    supply = get_power_supply(
        conf["mech"],
        conf["qdot_fuel"],
        conf["X_fuel"],
        conf["X_oxid"],
        conf["phase"],
    )

    D_fuel, D_oxid, e_wall = get_equivalent_diameters(
        conf["diam_fuel"],
        conf["diam_oxid_ext"],
        conf["diam_oxid_int"],
        conf["thick_wall"],
    )

    table_fuel = flow_workflow(
        name = "FUEL",
        mdot = supply.fuel_mass,
        rho0 = supply.fuel_normal_density,
        T    = 298.15,
        A    = np.pi * (D_fuel / 2)**2,
        save = save,
    )
    table_oxid = flow_workflow(
        name = "OXID",
        mdot = supply.oxidizer_mass,
        rho0 = supply.oxidizer_normal_density,
        T    = 273.15 + air_temp,
        A    = np.pi * (D_oxid / 2)**2,
        save = save,
    )

    return supply, table_fuel, table_oxid


def get_equivalent_diameters(
        diam_fuel: float,
        diam_oxid_ext: float,
        diam_oxid_int: float,
        thick_wall: float,
    ) -> tuple[float, float, float]:
    """ Compute equivalent diameters for fuel/oxidizer. """
    A_air_ext   = np.pi * (diam_oxid_ext / 2)**2
    A_air_int   = np.pi * (diam_oxid_int / 2)**2
    A_fuel_pipe = np.pi * (diam_fuel / 2 + thick_wall)**2
    A_fuel_tot  = 4 * A_fuel_pipe

    A_air = A_air_ext - A_air_int - A_fuel_tot
    A_air_one = A_air / 4

    e_wall = float(0.001 * thick_wall)
    D_fuel = float(0.001 * diam_fuel)
    D_oxid = float(0.001 * np.sqrt(4 * (A_air_one + A_fuel_pipe) / np.pi))

    return D_fuel, D_oxid, e_wall


def prepare_dimensions(save=False):
    """ Prepare geometry dimensions on the fly."""
    conf = load_yaml("initialize.yaml")

    D_fuel, D_oxid, e_wall = get_equivalent_diameters(
        conf["diam_fuel"],
        conf["diam_oxid_ext"],
        conf["diam_oxid_int"],
        conf["thick_wall"],
    )

    geometry_data = {
        "wedge_angle": 2.0,
        "len_inlet_ax": 3 * D_fuel,
        "len_inlet_an": 3 * D_fuel,
        "len_disperse": 1.00,
        "len_flue_out": 3.00,
        "dia_inlet_ax": D_fuel,
        "thk_inlet_wl": e_wall,
        "thk_inlet_an": D_oxid / 2,
        "thk_side_box": 2 * D_oxid,
    }

    if save:
        yaml = YAML()
        yaml.default_flow_style = False
        yaml.indent(mapping=2, sequence=4, offset=2)

        output_file = Path("domain.yaml")

        with output_file.open("w", encoding="utf-8") as f:
            yaml.dump(geometry_data, f)

    return geometry_data


def tabulate_dimensions():
    dims = prepare_dimensions()

    D_fuel     = 1000 * dims["dia_inlet_ax"]
    D_oxid_int = D_fuel + 2000 * dims["thk_inlet_wl"]
    D_oxid_ext = D_oxid_int + 2000 * dims["thk_inlet_an"]
    D_domain   = D_oxid_ext + 2000 * dims["thk_side_box"]
    L_domain   = dims["len_disperse"] + dims["len_flue_out"]

    return tabulate(
        [
            ("Fuel inlet diameter",           "mm", f"{D_fuel:.0f}"),
            ("Oxidizer inlet inner diameter", "mm", f"{D_oxid_int:.0f}"),
            ("Oxidizer inlet outer diameter", "mm", f"{D_oxid_ext:.0f}"),
            ("Domain diameter",               "mm", f"{D_domain:.0f}"),
            ("Domain length",                 "m",  L_domain),
        ],
        headers  = ["Quantity", "Unit", "Value"],
        tablefmt = "github"
    )


def mean_velocity(mdot, rho, T, A):
    """ Temperature corrected mean flow velocity. """
    return mdot / (rho * A)


def characteristic_length(A, f=0.07):
    """ Characteristic length for turbulence. """
    return f * np.sqrt(4.0 * A / np.pi)


def turbulent_kinetic_energy(U, I=0.05):
    """ Computes the turbulent kinetic energy. """
    return 1.5 * (U * I)**2


def turbulent_dissipation_rate(k, L, C_mu=0.09):
    """ Computes the turbulent dissipation rate. """
    return (C_mu**0.75) * (k**1.5) / L


def flow_workflow(name, mdot, rho0, T, A, I=0.05, save=False):
    rho = (273.15 / T) * rho0

    U = mean_velocity(mdot, rho, T, A)
    L = characteristic_length(A)
    k = turbulent_kinetic_energy(U, I)
    e = turbulent_dissipation_rate(k, L)

    if save:
        key_name = name.upper()
        _SETUP.set(f"{key_name}_T", T)
        _SETUP.set(f"{key_name}_U", U)
        _SETUP.set(f"{key_name}_RHO", rho)
        _SETUP.set(f"{key_name}_L", L)
        _SETUP.set(f"{key_name}_I", I)
        _SETUP.set(f"{key_name}_K", k)
        _SETUP.set(f"{key_name}_EPSILON", e)
        _SETUP.save()

    return tabulate(
        [
            ("Temperature", "K", T),
            ("Mean velocity", "m/s", U),
            ("Characteristic length", "mm", 1000*L),
            ("Turbulent kinetic energy", "m²/s²", k),
            ("Turbulent dissipation rate", "m²/s³", e),
        ],
        headers  = ["Quantity", "Unit", "Value"],
        tablefmt = "github"
    )


def read_data(
        case_file: str | Path = "case.foam",
        decomposed: bool = True,
        time_index: int = -1,
        verbose: bool = True,
    ) -> pv.DataObject | None:
    """ Read case data for post-processing. """
    case_file = Path(case_file)
    case_dir = case_file.parent

    if not case_file.exists():
        raise FileNotFoundError(f"No such case {case_file}")

    reader = pv.POpenFOAMReader(case_file)
    reader.enable_all_cell_arrays()
    reader.enable_all_point_arrays()

    if decomposed and list(case_dir.glob("processor*")):
        reader.case_type = "decomposed"

    try:
        case_time = reader.time_values[time_index]

        if verbose:
            print(f"Loading results for iter. {case_time}")

        reader.set_active_time_value(case_time)
        return reader.read()
    except AttributeError:
        print(f"Case was not ready to be read")
        return


def load_slice(*args, **kwargs) -> pv.PolyData | None:
    """ Standard slicing of mid-plane of domain. """
    if (mesh := kwargs.pop("mesh", None)) is None:
        if (mesh := read_data(*args, **kwargs)) is None:
            return

    return (
        mesh["internalMesh"]
            .slice(normal="z", origin=(0.0, 0.0, 0.0))
            .cell_data_to_point_data()
    )


def plot_slice(
        mesh,
        saveas: str | Path | None = None,
        off_screen: bool = True,
        window_size: tuple[int, int] = (900, 250),
        **kwargs
    ) -> None:
    """ Custom display of slice results. """
    # if off_screen and not saveas:
    #     print("Running headless without `saveas` has no effect, exit.")
    #     return

    kwargs.setdefault("cmap", "jet")
    kwargs.setdefault("scalar_bar_args", {
        "vertical": False,
        "height": 0.08,
        "width": 0.8,
        "position_x": 0.1,
        "title_font_size": 12,
        "label_font_size": 10,
    })

    pl = pv.Plotter(off_screen=off_screen)
    pl.add_mesh(mesh, **kwargs, show_edges=False)

    pl.camera_position = "xy"
    pl.window_size = window_size

    pl.enable_parallel_projection()
    pl.enable_2d_style()
    pl.enable_zoom_style()

    pl.add_ruler(
        pointa=[0.000, -0.02, 0.0],
        pointb=[4.001, -0.02, 0.0],
        title="Coordinate [m]"
    )
    region_bounds = [-0.1, 4.0, 0.0, 0.3, -0.01, 0.01]
    pl.view_xy(bounds=region_bounds)
    pl.zoom_camera(3.5)

    if saveas is not None:
        pl.show()
        pl.screenshot(saveas)

    return pl


def add_subplot(pl, mesh, **kwargs):
    kwargs.setdefault("show_edges", False)
    kwargs.setdefault("scalar_bar_args", {
        "vertical": False,
        "height": 0.1,
        "width": 0.8,
        "position_x": 0.1,
        "title_font_size": 12,
        "label_font_size": 10,
        "title": "",
    })

    xlim = kwargs.pop("xlim", (-0.1, 4.0))
    ylim = kwargs.pop("ylim", ( 0.0, 0.3))

    pl.add_mesh(mesh, show_scalar_bar=True, **kwargs)

    # TODO add this to the lessons learned.
    pl.scalar_bars.clear()

    pl.camera_position = "xy"
    pl.enable_parallel_projection()
    pl.enable_2d_style()
    pl.enable_zoom_style()

    pl.add_ruler(
        pointa=[0.000, -0.02, 0.0],
        pointb=[4.001, -0.02, 0.0],
        title="Coordinate [m]"
    )

    pl.view_xy(bounds=[*xlim, *ylim, -0.01, 0.01])
    pl.zoom_camera(2.85)


def add_suptitle(pl, title):
    pl.add_text(title, position="upper_left", font_size=12, color="black")
    # pl.add_text(title, position=(0.05, 0.96), font_size=12, color="black",
    #             viewport=True)


def plot_comparison(mesh1, mesh2, *, scalar, **kwargs):
    field_labels = {
        "T":   "Temperature [K]",
        "U":   "Mean velocity [m/s]",
        "O2":  "Oxygen mass fraction [-]",
        "CO":  "Carbon monoxide mass fraction [-]",
        "CO2": "Carbon dioxide mass fraction [-]",
        "H2O": "Water vapor mass fraction [-]",
        "a":   "Net absorption coefficient [1/m]"
    }

    label = field_labels.get(scalar, scalar)
    kwargs.pop("scalars", None)
    kwargs.setdefault("scalars", scalar)

    pl = pv.Plotter(shape=(2, 1), off_screen=True)
    pl.window_size = (900, 600)

    pl.subplot(0, 0)
    add_subplot(pl, mesh1, **kwargs)
    add_suptitle(pl, f"Cold air inlet - {label}")

    pl.subplot(1, 0)
    add_subplot(pl, mesh2, **kwargs)
    add_suptitle(pl, f"Pre-heated air inlet - {label}")

    pl.show()


@mj.plot(shape=(3, 1), size=(6, 8), sharex=True)
def plot_reports(*, plot, root="."):
    fig, ax = plot.subplots()
    post = mj.FoamPostProcessingLoader(root=root)

    df = post.load_report("outletT")
    x = df.iloc[:, 0].to_numpy()
    y = df.iloc[:, 1].to_numpy()

    ax[0].plot(x, y, color="k")
    ax[0].set_ylabel(f"Temperature [K]")

    df = post.load_report("outletCO")
    x = df.iloc[:, 0].to_numpy()
    y = df.iloc[:, 1].to_numpy() * 1e6

    ax[1].plot(x, y, color="k")
    ax[1].set_ylabel(f"Carbon monoxide [ppmm]")

    df = post.load_report("probe", select=r"**/T")
    x = df.iloc[:, 0].to_numpy()
    y = df.iloc[:, 1].to_numpy()

    ax[2].plot(x, y, color="k")
    ax[2].set_ylabel(f"Probe temperature [K]")

    for axis in ax:
        axis.grid(False)
        axis.set_xlabel("Iteration counter")


@mj.plot(shape=(3, 1), size=(12, 12), sharex=True)
def plot_fields(mesh, *, plot, **kwargs):
    x_points = [0.5, 1.0, 1.5, 2.0, 3.0, 3.5]
    resolution = kwargs.get("resolution", 200)

    def plot_field(idx, name, scale=None):
        scale = scale or (lambda x: x)

        for x_coord in x_points:
            line_data = mesh.sample_over_line(
                [x_coord, 0.00, 0.0],
                [x_coord, 0.40, 0.0],
                resolution=resolution
            )
            x = line_data["Distance"] * 100
            y = scale(line_data[name])
            ax[idx].plot(x, y, label=f"{x_coord:.2f} m")

    fig, ax = plot.subplots()

    plot_field(0, "T")
    ax[0].set_ylabel(f"Temperature [K]")

    plot_field(1, "U-normed")
    ax[1].set_ylabel(f"Velocity [m/s]")

    plot_field(2, "CO", scale=lambda x: 1e6 * x)
    ax[2].set_ylabel(f"Carbon monoxide [ppmm]")

    for axis in ax:
        axis.grid(False)
        axis.set_xlim(0.0, 40.0)
        axis.legend(loc=1, fontsize="small", ncol=2)
        axis.set_xlabel("Distance from axis [cm]")


def load_case(
        name,
        decomposed,
        show_plot = True,
        force_plot = False,
        **kwargs
    ):
    mesh = load_slice(
        case_file  = f"{name}/case.foam",
        decomposed = decomposed,
        **kwargs
    )

    if (decomposed and show_plot) or force_plot:
        plot_reports(root=name)

    return mesh


def handle_isosurfaces(pl, mesh, n_iso, scalar):
    if n_iso > 1:
        contour = mesh.contour(isosurfaces=n_iso, scalars=scalar)
        pl.add_mesh(contour, color="white", show_scalar_bar=False)

    pl.show()


def plot_temperature(mesh, n_iso=4):
    pl = plot_slice(mesh, scalars="T", cmap="hot")
    handle_isosurfaces(pl, mesh, n_iso, "T")


def plot_velocity(mesh):
    pl = plot_slice(mesh, scalars="U", cmap="jet")
    pl.show()


def plot_absorption(mesh, n_iso=5):
    pl = plot_slice(mesh, scalars="a", cmap="hot")
    handle_isosurfaces(pl, mesh, n_iso, "a")


def plot_oxygen(mesh, n_iso=5):
    pl = plot_slice(mesh, scalars="O2", cmap="coolwarm")
    handle_isosurfaces(pl, mesh, n_iso, "O2")


def plot_carbon_monoxide(mesh, n_iso=5):
    pl = plot_slice(mesh, scalars="CO", cmap="coolwarm")
    handle_isosurfaces(pl, mesh, n_iso, "CO")

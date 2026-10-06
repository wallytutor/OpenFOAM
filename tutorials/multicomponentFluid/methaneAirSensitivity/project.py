# -*- coding: utf-8 -*-

#region: to majordome
def import_headless_pyvista():
    """ Import PyVista in a headless environment. """
    import os

    # Mute C-level stderr (File Descriptor 2) immediately
    _devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(_devnull, 2)
    os.close(_devnull)

    # Force headless off-screen environment
    os.environ["VTK_DEFAULT_RENDER_WINDOW_HEADLESS"] = "1"
    os.environ["PYVISTA_OFF_SCREEN"] = "true"
    os.environ["PYVISTA_USE_IPYVISTA"] = "false"

    import pyvista as pv
    pv.set_jupyter_backend("static")

    return pv
#endregion

#region: imports
import cantera as ct
import majordome as mj
import numpy as np

from argparse import ArgumentParser
from pathlib import Path

from ruamel.yaml import YAML
from tabulate import tabulate

from IPython.display import Markdown, display

pv = import_headless_pyvista()
#endregion

#region: constants
# Molar composition of fuel [-].
X_FUEL = "CH4: 1"

# Molar composition of oxidizer [-].
X_OXID = "N2: 0.79, O2: 0.21"

# Mechanism and associated phase.
_MECHANISM = "chemkin/BFER_methane.yaml"
_PHASE     = "CH4_BFER_multi"

# Nominal fuel volume flow rate for 4 injectors [Nm³/h].
QDOT_FUEL = 500.0

# Base burner dimensions [m].
_DIAM_FUEL     = 0.027
_DIAM_OXID_EXT = 0.342
_DIAM_OXID_INT = 0.076
_THICK_WALL    = 0.0025

# Derived equivalent single-burner geometry and areas [m].
_A_AIR_EXT   = (np.pi / 4.0) * (_DIAM_OXID_EXT**2)
_A_AIR_INT   = (np.pi / 4.0) * (_DIAM_OXID_INT**2)
_A_FUEL_PIPE = (np.pi / 4.0) * ((_DIAM_FUEL + 2.0 * _THICK_WALL)**2)
_A_AIR_ONE   = (_A_AIR_EXT - _A_AIR_INT - 4.0 * _A_FUEL_PIPE) / 4.0

D_FUEL = _DIAM_FUEL
D_OXID = float(np.sqrt(4.0 * (_A_AIR_ONE + _A_FUEL_PIPE) / np.pi))
E_WALL = _THICK_WALL

AREA_FUEL = (np.pi / 4.0) * (D_FUEL**2)
AREA_OXIDIZER = (np.pi / 4.0) * (D_OXID**2)

# Domain dimensions [m].
_LEN_DISPERSE = 1.000
_LEN_FLUE_OUT = 3.000
_LEN_DOMAIN   = _LEN_DISPERSE + _LEN_FLUE_OUT

# Shared handle to create OpenFOAM initial conditions file.
_SETUP = mj.FoamDictFile(path="model/0/_shared")
#endregion

#region: problem setup
def get_power_supply(
        qdot_fuel: float = QDOT_FUEL,
        power: float | None = None,
        mech: str = _MECHANISM,
        phase: str = _PHASE,
        X_fuel: str | dict[str, float] = X_FUEL,
        X_oxid: str | dict[str, float] = X_OXID,
    ) -> mj.CombustionPowerSupply:
    """ Evaluate equivalent combustion power supply. """
    if power is None:
        normal_flow = mj.NormalFlowRate(mech, X=X_fuel, name=phase)
        mdot_fuel = normal_flow(qdot_fuel) / 4.0

        chon = mj.CombustionAtmosphereCHON(mech, phase=phase, basis="mole")
        lhv = chon.solution_heating_value(X_fuel, X_oxid)
        power = lhv * mdot_fuel * 1000.0

    supply = mj.CombustionPowerSupply(
        power       = power,
        equivalence = 1.0,
        fuel        = X_fuel,
        oxidizer    = X_oxid,
        mechanism   = mech,
        phase       = phase,
        basis       = "mole",
    )
    return supply


def flow_workflow(
        name: str,
        mdot: float,
        rho0: float,
        T: float,
        A: float,
        I: float = 0.05,
        f: float = 0.07,
        C_mu: float = 0.09,
        save: bool = False,
        show: bool = False,
    ) -> str:
    """ Evaluate turbulence parameters for a flow. """
    # Correct density for temperature.
    rho = (273.15 / T) * rho0

    # Temperature-corrected mean flow velocity.
    U = mdot / (rho * A)

    # Characteristic length for turbulence.
    L = f * np.sqrt(4.0 * A / np.pi)

    # Evaluate the turbulent kinetic energy.
    k = 1.5 * (U * I)**2

    # Evaluate the turbulent dissipation rate.
    e = (C_mu**0.75) * (k**1.5) / L

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

    table = tabulate(
        [
            ("Temperature", "K", T),
            ("Mean velocity", "m/s", U),
            ("Characteristic length", "mm", 1000.0 * L),
            ("Turbulent kinetic energy", "m²/s²", k),
            ("Turbulent dissipation rate", "m²/s³", e),
        ],
        headers  = ["Quantity", "Unit", "Value"],
        tablefmt = "github",
    )

    if show:
        display(Markdown(table))

    return table


def load_setup_shared(
        air_temp: float = 25.0,
        save: bool = False,
    ):
    """ Load simulation setup for the shared files. """
    supply = get_power_supply()

    table_fuel = flow_workflow(
        name = "FUEL",
        mdot = supply.fuel_mass,
        rho0 = supply.fuel_normal_density,
        T    = 298.15,
        A    = AREA_FUEL,
        save = save,
    )
    table_oxid = flow_workflow(
        name = "OXID",
        mdot = supply.oxidizer_mass,
        rho0 = supply.oxidizer_normal_density,
        T    = 273.15 + air_temp,
        A    = AREA_OXIDIZER,
        save = save,
    )

    return supply, table_fuel, table_oxid


def prepare_dimensions(save: bool = False) -> dict[str, float]:
    """ Prepare geometry dimensions on the fly. """
    geometry_data = {
        "wedge_angle": 2.0,
        "len_inlet_ax": 3 * D_FUEL,
        "len_inlet_an": 3 * D_FUEL,
        "len_disperse": _LEN_DISPERSE,
        "len_flue_out": _LEN_FLUE_OUT,
        "dia_inlet_ax": D_FUEL,
        "thk_inlet_wl": E_WALL,
        "thk_inlet_an": D_OXID / 2.0,
        "thk_side_box": 2.0 * D_OXID,
    }

    if save:
        yaml = YAML()
        yaml.default_flow_style = False
        yaml.indent(mapping=2, sequence=4, offset=2)

        output_file = Path("domain.yaml")

        with output_file.open("w", encoding="utf-8") as f:
            yaml.dump(geometry_data, f)

    return geometry_data
#endregion

#region: report only
def tabulate_dimensions(show: bool = False) -> str:
    """ Display dimensions of domain and inlets. """
    D_fuel     = 1000.0 * D_FUEL
    D_oxid_int = D_fuel + 2000.0 * E_WALL
    D_oxid_ext = D_oxid_int + 1000.0 * D_OXID
    D_domain   = D_oxid_ext + 4000.0 * D_OXID
    L_domain   = _LEN_DOMAIN

    table = tabulate(
        [
            ("Fuel inlet diameter",           "mm", f"{D_fuel:.0f}"),
            ("Oxidizer inlet inner diameter", "mm", f"{D_oxid_int:.0f}"),
            ("Oxidizer inlet outer diameter", "mm", f"{D_oxid_ext:.0f}"),
            ("Domain diameter",               "mm", f"{D_domain:.0f}"),
            ("Domain length",                 "m",  f"{L_domain:.0f}"),
        ],
        headers  = ["Quantity", "Unit", "Value"],
        tablefmt = "github",
    )

    if show:
        display(Markdown(table))

    return table
#endregion

#region: postprocessing
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
    except (AttributeError, IndexError):
        print("Case was not ready to be read")
        return None


def load_case(
        name: str,
        decomposed: bool = True,
        show_plot: bool = True,
        force_plot: bool = False,
        **kwargs
    ) -> pv.PolyData | None:
    """ Load OpenFOAM case slice and optionally plot convergence. """
    mesh = read_data(
        case_file  = f"{name}/case.foam",
        decomposed = decomposed,
        **kwargs
    )
    if mesh is None:
        return None

    mesh = (
        mesh["internalMesh"]
            .slice(normal="z", origin=(0.0, 0.0, 0.0))
            .cell_data_to_point_data()
    )

    if (decomposed and show_plot) or force_plot:
        plot_reports(root=name)

    return mesh


@mj.plot(shape=(3, 1), size=(6, 8), sharex=True)
def plot_reports(*, plot, root="."):
    fig, ax = plot.subplots()
    post = mj.FoamPostProcessingLoader(root=root)

    df = post.load_report("outletT")
    x = df.iloc[:, 0].to_numpy()
    y = df.iloc[:, 1].to_numpy()

    ax[0].plot(x, y, color="k")
    ax[0].set_ylabel("Temperature [K]")

    df = post.load_report("outletCO")
    x = df.iloc[:, 0].to_numpy()
    y = df.iloc[:, 1].to_numpy() * 1e6

    ax[1].plot(x, y, color="k")
    ax[1].set_ylabel("Carbon monoxide [ppmm]")

    df = post.load_report("probe", select=r"**/T")
    x = df.iloc[:, 0].to_numpy()
    y = df.iloc[:, 1].to_numpy()

    ax[2].plot(x, y, color="k")
    ax[2].set_ylabel("Probe temperature [K]")

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
    ax[0].set_ylabel("Temperature [K]")

    plot_field(1, "U-normed")
    ax[1].set_ylabel("Velocity [m/s]")

    plot_field(2, "CO", scale=lambda x: 1e6 * x)
    ax[2].set_ylabel("Carbon monoxide [ppmm]")

    for axis in ax:
        axis.grid(False)
        axis.set_xlim(0.0, 40.0)
        axis.legend(loc=1, fontsize="small", ncol=2)
        axis.set_xlabel("Distance from axis [cm]")


def plot_slice(
        mesh,
        window_size: tuple[int, int] = (900, 250),
        **kwargs
    ):
    """ Custom display of slice results. """
    kwargs.setdefault("cmap", "jet")
    kwargs.setdefault("scalar_bar_args", {
        "vertical": False,
        "height": 0.08,
        "width": 0.8,
        "position_x": 0.1,
        "title_font_size": 12,
        "label_font_size": 10,
    })

    pl = pv.Plotter(off_screen=True)
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
    return pl


def handle_isosurfaces(pl, mesh, n_iso, scalar):
    if n_iso > 1:
        contour = mesh.contour(isosurfaces=n_iso, scalars=scalar)
        pl.add_mesh(contour, color="white", show_scalar_bar=False)

    pl.show()


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
#endregion

#region: selected-plots
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
#endregion


def main():
    """ CLI workflow for case setup generation. """
    parser = ArgumentParser(
        description = "Prepare case conditions."
    )
    parser.add_argument(
        "--temperature",
        type    = float,
        default = 25.0,
        help    = "Air inlet temperature in degrees Celsius."
    )
    parser.add_argument(
        "--flow-rate",
        type    = float,
        default = QDOT_FUEL,
        help    = "Fuel volume flow rate in Nm³/h."
    )
    parser.add_argument(
        "--power",
        type    = float,
        default = None,
        help    = "Total power supply in kilo-watts."
    )
    args = parser.parse_args()

    if args.temperature < -273.15:
        print("Temperature must be above absolute zero.")
        return 1

    supply = get_power_supply(qdot_fuel=args.flow_rate, power=args.power)

    flow_workflow(
        name = "FUEL",
        mdot = supply.fuel_mass,
        rho0 = supply.fuel_normal_density,
        T    = 298.15,
        A    = AREA_FUEL,
        show = False,
        save = True,
    )

    flow_workflow(
        name = "OXID",
        mdot = supply.oxidizer_mass,
        rho0 = supply.oxidizer_normal_density,
        T    = args.temperature + 273.15,
        A    = AREA_OXIDIZER,
        show = False,
        save = True,
    )


if __name__ == "__main__":
    main()

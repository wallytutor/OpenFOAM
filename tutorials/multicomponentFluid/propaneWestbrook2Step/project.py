# -*- encoding: utf-8 -*-

#region: imports
import cantera as ct
import majordome as mj
import numpy as np

from argparse import ArgumentParser
from pathlib import Path

from numpy.typing import NDArray
from ruamel.yaml import YAML
from tabulate import tabulate

from IPython.display import Markdown, display

pv = mj.import_headless_pyvista()
#endregion

#region: constants
# Molar composition of fuel [-].
X_FUEL = "C3H8: 1"

# Molar composition of oxidizer [-].
X_OXID = "N2: 0.79, O2: 0.21"

# Mechanis and associated phase.
_MECHANISM = "chemkin/wb.yaml"
_PHASE     = "gas"

# System dimensions [m].
_DIAM_FUEL     = 0.027
_DIAM_OXID_INT = 0.032
_DIAM_OXID_EXT = 0.140
_DIAM_DOMAIN   = 1.000
_LEN_DOMAIN    = 4.000

# Shared handle to create initialConditions file.
_SETUP = mj.FoamDictFile(path="model/constant/initialConditions")

AREA_FUEL = (np.pi / 4.0) * (_DIAM_FUEL**2)
AREA_OXIDIZER = (np.pi / 4.0) * (_DIAM_OXID_EXT**2 - _DIAM_OXID_INT**2)
#endregion

#region: problem setup
def get_power_supply(
        power: float,
        mech: str = _MECHANISM,
        phase: str = _PHASE,
        X_fuel: str | dict[str, float] = X_FUEL,
        X_oxid: str | dict[str, float] = X_OXID,
    ) -> mj.CombustionPowerSupply:
    """ Evaluate equivalent combustion power supply. """
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

    if show:
        table = tabulate(
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
        display(Markdown(table))
#endregion

#region: report only
def tabulate_dimensions():
    """ Display dimensions of domain and inlets. """
    D_fuel     = 1000 * _DIAM_FUEL
    D_oxid_int = 1000 * _DIAM_OXID_INT
    D_oxid_ext = 1000 * _DIAM_OXID_EXT
    D_domain   = 1000 * _DIAM_DOMAIN
    L_domain   = _LEN_DOMAIN

    table = tabulate(
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
    display(Markdown(table))
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
    except AttributeError:
        print(f"Case was not ready to be read")
        return


def load_case(
        name: str,
        decomposed: bool = True,
        show_plot: bool = True,
        force_plot: bool = False,
        **kwargs
    ):
    mesh = read_data(
        case_file  = f"{name}/case.foam",
        decomposed = decomposed,
        **kwargs
    )

    mesh = mesh["internalMesh"] \
            .slice(normal="z", origin=(0.0, 0.0, 0.0)) \
            .cell_data_to_point_data()

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
    ax[0].set_ylabel(f"Outlet temperature [K]")

    df = post.load_report("outletCO")
    x = df.iloc[:, 0].to_numpy()
    y = df.iloc[:, 1].to_numpy() * 1e6

    ax[1].plot(x, y, color="k")
    ax[1].set_ylabel(f"Outlet CO [ppmm]")

    df = post.load_report("probe", select=r"**/T")
    x = df.iloc[:, 0].to_numpy()

    for col_name in df.columns[1:]:
        y = df[col_name].to_numpy()
        ax[2].plot(x, y, label=col_name)

    ax[2].set_ylabel(f"Probe temperature [K]")
    ax[2].legend(loc=2, fontsize="x-small")

    for axis in ax:
        axis.grid(False)
        axis.set_xlabel("Iteration counter")


@mj.plot(shape=(3, 1), size=(12, 12), sharex=True)
def plot_fields(mesh, *, plot, **kwargs):
    x_points = [0.5, 1.0, 1.5, 2.0, 3.0, 3.5]
    resolution = kwargs.get("resolution", 200)
    x_max = _DIAM_DOMAIN / 2

    def plot_field(idx, name, scale=None):
        scale = scale or (lambda x: x)

        for x_coord in x_points:
            line_data = mesh.sample_over_line(
                [x_coord, 0.00, 0.0],
                [x_coord, 0.999 * x_max, 0.0],
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
        axis.set_xlim(0.0, 100 * x_max)
        axis.legend(loc=1, fontsize="small", ncol=2)
        axis.set_xlabel("Distance from axis [cm]")


def evaluate_arrows(mesh, nx: int, ny: int, factor: float, padding: float):
    # (xmin, xmax, ymin, ymax, zmin, zmax)
    bounds = mesh.bounds

    # Define a small inward margin
    margin_x = (bounds[1] - bounds[0]) * padding
    margin_y = (bounds[3] - bounds[2]) * padding

    # Create a slightly shrunken uniform grid
    grid = pv.ImageData(
        dimensions=(nx, ny, 1),
        spacing=(
            (bounds[1] - bounds[0] - 2 * margin_x) / (nx - 1),
            (bounds[3] - bounds[2] - 2 * margin_y) / (ny - 1),
            1.0
        ),
        origin=(bounds[0] + margin_x, bounds[2] + margin_y, bounds[4])
    )

    # Resample and generate fixed-size arrows
    resampled_grid = grid.sample(mesh)
    resampled_grid.set_active_vectors("U")
    arrows = resampled_grid.glyph(scale=False, orient="U", factor=factor)

    # mesh.set_active_vectors("U")
    # arrows = mesh.glyph(
    #     scale     = False,
    #     orient    = "U",
    #     factor    = vector_factor,
    #     tolerance = vector_tolerance
    # )

    return arrows


def plot_slice(
        mesh,
        window_size: tuple[int, int] = (900, 250),
        show_vectors: bool = True,
        factor: float = 0.04,
        padding: float = 0.025,
        nx: int = 45,
        ny: int = 10,
        **kwargs
    ) -> None:
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

    if show_vectors:
        arrows = evaluate_arrows(mesh, nx, ny, factor, padding)
        pl.add_mesh(arrows, color="black")

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


def plot_propane(mesh, n_iso=5):
    pl = plot_slice(mesh, scalars="C3H8", cmap="coolwarm", clim=[0.0, 0.01])
    handle_isosurfaces(pl, mesh, n_iso, "C3H8")


def plot_oxygen(mesh, n_iso=5):
    pl = plot_slice(mesh, scalars="O2", cmap="coolwarm")
    handle_isosurfaces(pl, mesh, n_iso, "O2")


def plot_carbon_monoxide(mesh, n_iso=5):
    pl = plot_slice(mesh, scalars="CO", cmap="coolwarm")
    handle_isosurfaces(pl, mesh, n_iso, "CO")


def plot_reciprocal_time_step(mesh):
    pl = plot_slice(mesh, scalars="rDeltaT", cmap="jet", show_vectors=False)
    pl.show()
#endregion

def main():
    """ CLI workflow for case setup generation. """
    parser = ArgumentParser(
        description = "Prepare case conditions."
    )
    parser.add_argument(
        "--power",
        type     = float,
        help     = "Total power supply in kilo-watts."
    )
    parser.add_argument(
        "--temperature",
        type     = float,
        help     = "Air inlet temperature in degrees Celsius."
    )
    args = parser.parse_args()

    if args.power is None or args.power <= 0.0:
        print("Power must be a positive value.")
        return 1

    if args.temperature is None or args.temperature <= 0.0:
        print("Temperature must be a positive value.")
        return 1

    supply = get_power_supply(args.power)

    flow_workflow(
        name = "fuel",
        mdot = supply.fuel_mass,
        rho0 = supply.fuel_normal_density,
        T    = 298.15,
        A    = AREA_FUEL,
        show = False,
        save = True,
    )

    flow_workflow(
        name = "oxid",
        mdot = supply.oxidizer_mass,
        rho0 = supply.oxidizer_normal_density,
        T    = args.temperature + 273.15,
        A    = AREA_OXIDIZER,
        show = False,
        save = True,
    )


if __name__ == "__main__":
    main()

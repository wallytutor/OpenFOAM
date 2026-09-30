# -*- coding: utf-8 -*-

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

_SETUP = mj.FoamDictFile(path="0.orig/_shared")


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


def get_equivalent_diameters(
        diam_fuel: float,
        diam_oxid_ext: float,
        diam_oxid_int: float,
        thick_wall: float,
    ) -> tuple[float, float]:
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

    yaml = YAML()
    yaml.default_flow_style = False
    yaml.indent(mapping=2, sequence=4, offset=2)

    output_file = Path("domain.yaml")

    with output_file.open("w", encoding="utf-8") as f:
        yaml.dump(geometry_data, f)

    return D_fuel, D_oxid


def fit_transport(
        mech: str,
        phase: str,
        T_space: tuple[float, float, int] = (400, 2300, 100)
    ) -> mj.SutherlandFitting:
    """ Fit Sutherland coefficients for a mechanism and phase. """
    sutherland = mj.SutherlandFitting(mech, name=phase)
    sutherland.fit(np.linspace(*T_space))

    transport = Path(mech).parent / "transport"
    transport.unlink(missing_ok=True)

    _ = sutherland.as_openfoam_dict(transport)

    # Use the following for inspection:
    # figs = Path(mech).parent / "media"
    # figs.mkdir(exist_ok=True)

    # for s in coef["species"]:
    #     p = sutherland.plot_species(s)
    #     p.savefig(figs / f"{s}.png")

    return sutherland


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


def flow_workflow(name, mdot, rho0, T, A, I=0.05):
    rho = (273.15 / T) * rho0

    U = mean_velocity(mdot, rho, T, A)
    L = characteristic_length(A)
    k = turbulent_kinetic_energy(U, I)
    e = turbulent_dissipation_rate(k, L)

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
        time_index: int = -1
    ) -> pv.DataObject | None:
    """ Read case data for post-processing. """
    if not Path(case_file).exists():
        raise FileNotFoundError(f"No such case {case_file}")

    reader = pv.POpenFOAMReader(case_file)
    reader.enable_all_cell_arrays()
    reader.enable_all_point_arrays()

    if decomposed and list(Path(".").glob("processor*")):
        reader.case_type = "decomposed"

    try:
        case_time = reader.time_values[time_index]
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
        off_screen: bool = False,
        window_size: tuple[int, int] = (900, 250),
        **kwargs
    ) -> None:
    """ Custom display of slice results. """
    if off_screen and not saveas:
        print("Running headless without `saveas` has no effect, exit.")
        return

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


@mj.plot(shape=(3, 1), size=(8, 9), sharex=True)
def plot_reports(*, plot, root="."):
    fig, ax = plot.subplots()
    post = mj.FoamPostProcessingLoader(root=root)

    df = post.load_report("outletT")
    x = df.iloc[:, 0].to_numpy()
    y = df.iloc[:, 1].to_numpy()

    ax[0].plot(x, y, color="k")
    ax[0].set_xlabel("Iteration counter")
    ax[0].set_ylabel(f"Temperature [K]")

    df = post.load_report("outletCO")
    x = df.iloc[:, 0].to_numpy()
    y = df.iloc[:, 1].to_numpy()

    ax[1].plot(x, y, color="k")
    ax[1].set_xlabel("Iteration counter")
    ax[1].set_ylabel(f"Carbon monoxide [-]")

    df = post.load_report("probe", select=r"**/T")
    x = df.iloc[:, 0].to_numpy()
    y = df.iloc[:, 1].to_numpy()

    ax[2].plot(x, y, color="k")
    ax[2].set_xlabel("Iteration counter")
    ax[2].set_ylabel(f"Probe temperature [K]")

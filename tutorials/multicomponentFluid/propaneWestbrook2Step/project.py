# -*- encoding: utf-8 -*-

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

import cantera as ct
import majordome as mj
import numpy as np

from argparse import ArgumentParser
from pathlib import Path

from numpy.typing import NDArray
from ruamel.yaml import YAML
from tabulate import tabulate

from IPython.display import Markdown, display

pv = import_headless_pyvista()

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
_DIAM_OXID_EXT = 0.202
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
        T    = args.air_temperature + 273.15,
        A    = AREA_OXIDIZER,
        show = False,
        save = True,
    )


if __name__ == "__main__":
    main()

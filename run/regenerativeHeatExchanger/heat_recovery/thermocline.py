# -*- coding: utf-8 -*-
import cantera as ct
import numpy as np
import sympy as sp
from majordome.utilities import MajordomePlot
from IPython.display import Markdown
from itertools import cycle
from matplotlib import cm
from ruamel.yaml import YAML


class SymbolicTracker:
    """ Utility class to track symbolic variables. """
    SYMBOLS = []

    @classmethod
    def get(cls, name):
        for sym in cls.SYMBOLS:
            if sym.name == name:
                return sym

        raise ValueError(f"Symbol with name '{name}' not found.")

    @classmethod
    def clear(cls):
        cls.SYMBOLS.clear()

    @classmethod
    def sym(cls, name, **kwargs) -> sp.Symbol:
        if name in [s.name for s in cls.SYMBOLS]:
            raise ValueError(f"Symbol '{name}' already exists.")

        symbol = sp.symbols(name, **kwargs)
        cls.SYMBOLS.append(symbol)
        return symbol

    @classmethod
    def lambdify(cls, expr, symbols=None, use_all=True):
        # Only use all symbols if none are provided
        if not symbols and use_all:
            symbols = cls.SYMBOLS

        if symbols is None:
            symbols = list(expr.free_symbols)
            symbols = sorted(symbols, key=lambda x: x.name)

        return sp.lambdify(symbols, expr, modules="numpy")


class SolidMaterial:
    """ Solid material properties. """
    def __init__(self):
        self._r = SymbolicTracker.sym("rho_s")
        self._c = SymbolicTracker.sym("c_ps")
        self._k = SymbolicTracker.sym("k_s")

    @property
    def density(self) -> sp.Symbol:
        """ Density of solid material [kg/m^3]. """
        return self._r

    @property
    def specific_heat_capacity(self) -> sp.Symbol:
        """ Specific heat capacity of solid material [J/kg-K]. """
        return self._c

    @property
    def thermal_conductivity(self) -> sp.Symbol:
        """ Thermal conductivity of solid material [W/m-K]. """
        return self._k

    def specific_energy(self, delta_T) -> sp.Symbol:
        """ Specific energy of solid material [J/m³]. """
        return self._r * self._c * delta_T  # type: ignore


class OperatingFluid:
    """ Operating fluid properties. """
    def __init__(self):
        self._h_c = SymbolicTracker.sym("h_g_c", positive=True)
        self._h_h = SymbolicTracker.sym("h_g_h", positive=True)
        self._r   = SymbolicTracker.sym("rho_g", positive=True)
        self._mu  = SymbolicTracker.sym("mu_g", positive=True)

    @property
    def density_inlet(self) -> sp.Symbol:
        """ Density of operating fluid at inlet [kg/m³]. """
        return self._r

    @property
    def enthalpy_cold(self) -> sp.Symbol:
        """ Enthalpy of cold fluid [J/kg]. """
        return self._h_c

    @property
    def enthalpy_hot(self) -> sp.Symbol:
        """ Enthalpy of hot fluid [J/kg]. """
        return self._h_h

    @property
    def specific_energy_inlet(self) -> sp.Symbol:
        """ Specific energy of operating fluid at inlet [J/m³]. """
        return self._r * (self._h_h - self._h_c)  # type: ignore

    @property
    def dynamic_viscosity(self):
        """ Dynamic viscosity of operating fluid [Pa-s]. """
        return self._mu


class HexagonalCell:
    """ Hexagonal cell *building block* of the thermocline model. """
    def __init__(self):
        self._D_h = SymbolicTracker.sym("D_h", positive=True)
        self._m_h = SymbolicTracker.sym("m_h", positive=True)

        # Half the distance between holes [m]
        l = (1 + self._m_h) * self._D_h / 2  # type: ignore

        # Height of hexagonal cell [m]
        h_h = sp.sqrt(l**2 - (l/2)**2)

        # Area of triangle formed by center of holes and symmetry [m^2]
        a_t = h_h * l / 2

        # Area of hexagonal cell [m^2]
        self._a_c = 6 * a_t

        # Area of holes in hexagonal cell [m^2]
        self._a_h = sp.pi * self._D_h**2 / 4  # type: ignore

        # Area of solid material in hexagonal cell [m^2]
        self._a_s = self._a_c - self._a_h  # type: ignore

    @property
    def solid_cross_section(self):
        """ Cross-section of solid material in hexagonal cell [m²]. """
        return self._a_s

    @property
    def fluid_cross_section(self):
        """ Cross-section of holes in hexagonal cell [m²]. """
        return self._a_h

    @property
    def surface_to_volume_ratio(self):
        """ Ratio of exchange area per unit volume [1/m]. """
        return sp.pi * self._D_h / self._a_s

    def get_surface_to_volume_ratio(self):
        expr = self.surface_to_volume_ratio
        return sp.lambdify((self._D_h, self._m_h), expr, "numpy")


class ThermoclineDescription:
    """ Thermocline model based on target operating conditions. """
    def __init__(self):
        self._H_t = SymbolicTracker.sym("H_t", positive=True)
        self._P_t = SymbolicTracker.sym("P_t", positive=True)
        self._h_t = SymbolicTracker.sym("h_t", positive=True)
        self._T_h = SymbolicTracker.sym("T_h", positive=True)
        self._T_c = SymbolicTracker.sym("T_c", positive=True)

    @property
    def height(self):
        """ Height of thermocline structure [m]. """
        return self._h_t

    @property
    def power(self):
        """ Target nominal operating power [W]. """
        return self._P_t

    @property
    def stored_energy(self):
        """ Target total stored energy [J]. """
        return self._H_t

    @property
    def temperature_hot(self):
        """ Hot working temperature [K]. """
        return self._T_h

    @property
    def temperature_cold(self):
        """ Cold working temperature [K]. """
        return self._T_c

    @property
    def temperature_difference(self):
        """ Temperature difference between hot and cold operation [K]. """
        return self._T_h - self._T_c  # type: ignore

    @property
    def specific_energy(self):
        """ Specific energy of thermocline model [J/m]. """
        return self._H_t / self._h_t  # type: ignore


class PressureDropCalculator:
    def __init__(self, L, D):
        self._L = L
        self._D = D

    def darcy_weisbach_pressure_drop(self, rho, v, f=0.05):
        """ Calculate pressure drop using the Darcy-Weisbach equation."""
        if callable(f):
            f = f(rho, v, self._D)

        return f * (self._L / self._D) * (rho * v**2 / 2)


class ThermoclineModel:
    def __init__(self, fname):
        self._load_symbols()
        self._load_data(fname)

    def _load_symbols(self):
        SymbolicTracker.clear()

        ###
        # Define models
        ###

        self.thermocline = thermocline = ThermoclineDescription()
        self.cell_block  = cell_block  = HexagonalCell()
        self.solid       = solid       = SolidMaterial()
        self.fluid       = fluid       = OperatingFluid()

        ###
        # General calculations
        ###

        self.S_V = cell_block.surface_to_volume_ratio
        self.A_t = total_cross_sectional_area(solid, thermocline)
        self.N_c = cells_per_cross_section(self.A_t, cell_block)
        self.P_c = power_per_cell(thermocline, self.N_c)
        self.U_g = fluid_injection_velocity(self.P_c, fluid, cell_block)
        self.Re  = reynolds_number(self.U_g, cell_block, fluid)

        ###
        # Organize symbols
        ###

        variables = [
            # ThermoclineDescription
            SymbolicTracker.get("H_t"),
            SymbolicTracker.get("h_t"),
            SymbolicTracker.get("T_h"),
            SymbolicTracker.get("T_c"),
            # HexagonalCell
            SymbolicTracker.get("D_h"),
            SymbolicTracker.get("m_h"),
            # SolidMaterial
            SymbolicTracker.get("rho_s"),
            SymbolicTracker.get("c_ps"),
            SymbolicTracker.get("k_s"),
            # OperatingFluid
            SymbolicTracker.get("h_g_c"),
            SymbolicTracker.get("h_g_h"),
            # Charging/discharging parameters
            SymbolicTracker.get("P_t"),
            SymbolicTracker.get("rho_g"),
            SymbolicTracker.get("mu_g"),
        ]

        ###
        # Lambdify functions
        ###

        def lambdify(expr):
            return SymbolicTracker.lambdify(expr, symbols=variables)

        self.fn_S_V = lambdify(self.S_V)
        self.fn_A_t = lambdify(self.A_t)
        self.fn_N_c = lambdify(self.N_c)
        self.fn_P_c = lambdify(self.P_c)
        self.fn_U_g = lambdify(self.U_g)
        self.fn_Re  = lambdify(self.Re)

    def _load_data(self, fname):
        MW = 1e6
        MWh = 3600 * MW

        ###
        # Load numerical values from YAML file
        ###

        yaml = YAML()
        parameters = yaml.load(open(fname))

        self.num_h_t   = parameters["h_t"]
        self.num_H_t   = parameters["H_t"] * MWh
        self.num_P_c   = parameters["P_c"] * MW
        self.num_P_d   = parameters["P_d"] * MW
        self.num_T_h   = parameters["T_h"]
        self.num_T_c   = parameters["T_c"]
        self.num_m_h   = parameters["m_h"]
        self.num_D_h   = parameters["D_h"]
        self.num_rho_s = parameters["rho_s"]
        self.num_c_ps  = parameters["c_ps"]
        self.num_k_s   = parameters["k_s"]

        args_c = air_properties(self.num_T_c)
        args_h = air_properties(self.num_T_h)

        self.num_rho_g_c, self.num_h_g_c, self.num_mu_g_c = args_c
        self.num_rho_g_h, self.num_h_g_h, self.num_mu_g_h = args_h

        ###
        # Organize argument groups
        ###

        args_all = [
            # ThermoclineDescription
            self.num_H_t,
            self.num_h_t,
            self.num_T_h,
            self.num_T_c,
            # HexagonalCell
            self.num_D_h,
            self.num_m_h,
            # SolidMaterial
            self.num_rho_s,
            self.num_c_ps,
            self.num_k_s,
            # OperatingFluid
            self.num_h_g_c,
            self.num_h_g_h
        ]

        self.args_charging = (
            *args_all, self.num_P_c, self.num_rho_g_h, self.num_mu_g_h)

        self.args_discharging = (
            *args_all, self.num_P_d, self.num_rho_g_c, self.num_mu_g_c)

        # Assume 0.005 as a representative value for the relative roughness of
        # refractory bricks, which have typical relative roughness in the range
        # from 0.001 to 0.01. The friction factor can be estimated using the
        # Swamee-Jain approximation, a simplified version of the implicit
        # Colebrook-White equation given by:
        # f = \dfrac{0.25}{\left[
        # \log_{10}\left(\dfrac{\epsilon}{3.7D} + \dfrac{5.74}{Re^{0.9}}\right)
        # \right]^2}
        # XXX check use, too high for transitional flow!
        # f = tc.swamee_jain_friction_factor(fn_Re(*args), 0.008, num_D_h)

        # Surface-to-volume ratio of the hexagonal cell:
        self.num_S_V = self.fn_S_V(*self.args_charging)

        # Diameter of the system:
        A = self.fn_A_t(*self.args_charging)
        self.num_D_sys = np.sqrt(4 * A / np.pi)

    def tabulate(self):
        rows = [
            ("Parameter",  "Value", "Units"),
            ("---------",  "-----", "-----"),

            ("Stored energy in thermocline",
             f"{self.num_H_t/3.6e9:.1f}", "MWh"),

            ("Target nominal operating power (charge)",
             f"{self.num_P_c/1e6:.1f}", "MW"),

            ("Target nominal operating power (discharge)",
             f"{self.num_P_d/1e6:.1f}", "MW"),

            ("Height of thermocline",
             f"{self.num_h_t:.1f}", "m"),

            ("Diameter of system",
             f"{self.num_D_sys:.1f}", "m"),

            ("Hot working temperature",
             f"{self.num_T_h:.1f}", "K"),

            ("Cold working temperature",
             f"{self.num_T_c:.1f}", "K"),

            ("Diameter of holes in hexagonal cell",
             f"{self.num_D_h*100:.1f}", "cm"),

            ("Wall thickness scaling factor",
             f"{self.num_m_h:.1f}", "-"),

            ("Density of solid material",
             f"{self.num_rho_s:.0f}", "kg/m³"),

            ("Specific heat capacity of solid material",
             f"{self.num_c_ps:.1f}", "J/kg-K"),

            ("Thermal conductivity of solid material",
             f"{self.num_k_s:.1f}", "W/m-K"),

            ("Enthalpy of cold fluid",
             f"{self.num_h_g_c/1e6:.3f}", "MJ/kg"),

            ("Enthalpy of hot fluid",
             f"{self.num_h_g_h/1e6:.3f}", "MJ/kg"),

            ("Density of operating fluid (cold)",
             f"{self.num_rho_g_c:.2f}", "kg/m³"),

            ("Density of operating fluid (hot)",
             f"{self.num_rho_g_h:.2f}", "kg/m³"),

            ("Dynamic viscosity of operating fluid",
             f"{self.num_mu_g_h*1e6:.2f}", "uPa-s"),

            ("Surface-to-volume ratio of hexagonal cell",
             f"{self.num_S_V:.3f}", "1/m")
        ]

        def format_row(name, value, units):
            return f"| {name} | {value} | {units} |\n"

        rows = "".join(format_row(*args) for args in rows)
        return Markdown(rows)


def total_cross_sectional_area(solid, thermocline):
    """ Total cross-sectional area of solid material in thermocline [m²]. """
    hv = solid.specific_energy(thermocline.temperature_difference)
    return thermocline.specific_energy / hv  # type: ignore


def cells_per_cross_section(A_t, cell_block):
    """ Number of hexagonal cells per cross section of thermocline. """
    return A_t / cell_block.solid_cross_section  # type: ignore


def fluid_injection_velocity(P_n, fluid, cell_block):
    """ Fluid velocity through holes in hexagonal cell [m/s]. """
    return P_n / (fluid.specific_energy_inlet * cell_block.fluid_cross_section)  # type: ignore


def reynolds_number(U_g, cell_block, fluid):
    """ Reynolds number for fluid flow through holes in hexagonal cell. """
    D = sp.sqrt(4 * cell_block.fluid_cross_section / sp.pi)  # type: ignore
    return fluid.density_inlet * U_g * D / fluid.dynamic_viscosity  # type: ignore


def swamee_jain_friction_factor(Re, epsilon, D):
    """ Calculate friction factor using the Swamee-Jain equation. """
    return 0.25 / (np.log10(epsilon / (3.7 * D) + 5.74 / Re**0.9))**2


def power_per_cell(thermocline, N):
    """ Power per hexagonal cell [W]. """
    return thermocline.power / N  # type: ignore


def air_properties(T):
    air = ct.Solution("airish.yaml")
    air.X = "O2:0.21, N2:0.79"
    air.TP = T, ct.one_atm
    return air.density_mass, air.enthalpy_mass, air.viscosity


@MajordomePlot.new(size=(5, 4), xlabel="Diameter of holes (cm)",
                   ylabel="Surface-to-volume ratio (1/m)")
def plot_surface_to_volume_ratio(cell_block, *, plot):
    func = cell_block.get_surface_to_volume_ratio()
    D_samples = np.linspace(0.010, 0.050, 100)
    m_samples = np.arange(0.5, 2.0+0.001, 0.5)

    _, ax = plot.subplots()
    ax = ax[0]

    color_map = cm.get_cmap("tab10")
    color_cycle = cycle(color_map.colors) # type: ignore

    for m in m_samples:
        ratio_sv = func(D_samples, m)
        color = next(color_cycle)
        ax.plot(100*D_samples, ratio_sv, color=color, label=f"m={m:.1f}")

    ax.legend(loc=1, fontsize="small")
    return plot


@MajordomePlot.new(size=(5, 4), xlabel="Diameter of holes (cm)",
                   ylabel="Pressure drop (mbar)")
def plot_pressure_drop(func, *, plot):
    D_samples = np.linspace(0.010, 0.050, 100)
    m_samples = np.arange(0.5, 2.0+0.001, 0.5)

    _, ax = plot.subplots()
    ax = ax[0]

    color_map = cm.get_cmap("tab10")
    color_cycle = cycle(color_map.colors) # type: ignore

    for m in m_samples:
        val = func(D_samples, m) / 100
        color = next(color_cycle)
        ax.plot(100*D_samples, val, color=color, label=f"m={m:.1f}")

    ax.legend(loc=1, fontsize="small")
    return plot

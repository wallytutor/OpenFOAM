# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from majordome.utilities import plot_xy


def conductivity(T):
    """ Reference thermal conductivity of material A. """
    return 1.9847 - T * (2.531209E-03 - 1.43E-06 * T)


def load_conductivity_data(fname="k_imperial.dat"):
    """ Load the conductivity data from the given file. """
    df = pd.read_csv(fname, sep=r"\s+", header=None, comment="#")
    # df.columns = ["T", "k_A", "k_B"]

    df[0] += 273.15
    df[1] *= 0.144131  # Btu/(hr·ft·°F) to W/(m·K)
    df[2] *= 0.144131  # Btu/(hr·ft·°F) to W/(m·K)

    return df


def plot_conductivity():
    """ Plot the conductivity data from the given dataframe. """
    df = load_conductivity_data()
    T = np.linspace(df[0].min(), df[0].max(), 100)
    k = conductivity(T)

    p = plot_xy()
    p.axes[0].plot(df[0].to_numpy(), df[1].to_numpy(), "ko", label="RS")
    p.axes[0].plot(df[0].to_numpy(), df[2].to_numpy(), "ro", label="AS")
    p.axes[0].plot(T, k, "k", label="_none_")
    p.axes[0].set_xlabel("Temperature [K]")
    p.axes[0].set_ylabel("Thermal conductivity [W/(m·K)]")
    p.axes[0].legend(loc=3, fontsize="small")
    p.savefig("../conductivity.png")


def print_table(x, y):
    """ Print a table of the given x and y values. """
    for xk, yk in zip(x, y):
        print(f"({xk:6.1f} {yk:6.3f})")


def tabulate_specific_heat_capacity(T):
    """ Tabulate the specific heat capacity of the material. """
    df = pd.read_csv("cp.dat", sep=r"\s+", header=None)
    cp = np.interp(T, df[0].to_numpy(),  df[1].to_numpy())

    print("\nSpecific heat capacity [J/(kg·K)]")
    print_table(T, cp)


def tabulate_conductivity(T, n):
    """ Tabulate the thermal conductivity of the material. """
    df = load_conductivity_data()
    kc = np.interp(T, df[0].to_numpy(),  df[n].to_numpy())

    print(f"\nThermal conductivity {n} [W/(m·K)]")
    print_table(T, kc)


if __name__ == "__main__":
    T_new = np.arange(200, 1301, 100)

    tabulate_specific_heat_capacity(T_new)
    tabulate_conductivity(T_new, n=1)
    tabulate_conductivity(T_new, n=2)
    plot_conductivity()

    # XXX: only rhoConst is supported, not worth printing...
    # df = pd.read_csv("rho.dat", sep=r"\s+", header=None)
    # rho_new = np.interp(T_new, df[0].to_numpy(),  df[1].to_numpy())

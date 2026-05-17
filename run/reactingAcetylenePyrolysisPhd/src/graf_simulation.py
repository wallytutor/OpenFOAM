# -*- coding: utf-8 -*-
from pathlib import Path
from tempfile import NamedTemporaryFile
from subprocess import run
from textwrap import dedent
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


def add_directory():
    """ Add path at repository root to Cantera path. """
    here = Path(__file__).resolve().parent
    ctdir =  here / "../../../../references/kinetics/graf_2007_global"
    ct.add_directory(ctdir)


class GrafSimulation:
    """ Simulate single reactor with Graf's mechanism. """
    def __init__(self, *, T, P, f, verbose=True):
        add_directory()

        gas = ct.Solution("graf-2007.yaml")
        gas.TPX = T, P, self._acetylene_mixture(f)

        reactor = ct.IdealGasConstPressureReactor(gas)
        reactor.energy_enabled = False

        simulation = ct.ReactorNet([reactor])
        simulation.verbose = verbose

        states = ct.SolutionArray(gas, extra=["t"])

        self._reactor = reactor
        self._simulation = simulation
        self._states = states

        self._feed_state()

    def _acetylene_mixture(self, f):
        """ Simulate single reactor with Graf's mechanism. """
        X = {"C2H2": 0.98 * f, "CH4": 0.002 * f}
        X["N2"] = 1.0 - sum(X.values())
        return X

    def _feed_state(self):
        """ Feed solution state to solution array. """
        state = self._reactor.thermo.state
        self._states.append(state, t=self._simulation.time)

    def run(self, *, tout, nsteps, ttol=1.0e-07):
        """ Simulate system until `tout` seconds in `nsteps` steps. """
        times = np.linspace(0, tout, nsteps)

        for tend in times[1:]:
            if self._simulation.time > tend:
                raise RuntimeError("Simulation time greater than exit time...")

            self._simulation.advance(tend)
            self._feed_state()

            if abs(tend - self._simulation.time) > ttol:
                raise RuntimeError("Unable to reach exit time during advance...")

        return self._states

    def plot(self):
        """ Displays model solution over time. """
        # Same indexing as Octave
        t  = self._states.t
        y1 = 100 * self._states.C2H2
        y2 = 100 * self._states.H2
        y3 = 100 * self._states.C2H4
        y4 = 100 * self._states.CH4
        y5 = 100 * self._states.C4H4
        y6 = 100 * self._states.C6H6
        y7 = 100 * self._states.C

        plt.close("all")
        fig, ax = plt.subplots(2, 2, figsize=(10, 8))
        ax = np.ravel(ax)

        def set_subplot(ax, ymin=0.0):
            ax.grid(linestyle=":")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Mass percentage")
            ax.legend(loc="best")
            ax.set_xlim(0, t[-1])
            ax.set_ylim(ymin, None)

        ax[0].plot(t, y1, c="b", label="$C_2H_2$")
        ax[1].plot(t, y2, c="b", label="$H_2$")
        ax[1].plot(t, y7, c="r", label="$C_s$")
        ax[2].plot(t, y3, c="b", label="$C_2H_4$")
        ax[2].plot(t, y4, c="r", label="$CH_4$")
        ax[3].plot(t, y5, c="b", label="$C_4H_4$")
        ax[3].plot(t, y6, c="r", label="$C_6H_6$")

        set_subplot(ax[0])
        set_subplot(ax[1])
        set_subplot(ax[2])
        set_subplot(ax[3])

        fig.tight_layout()
        return fig, ax


def make_octave_program(states, saveas="graf_plot_octave"):
    """ Generate an Octave program with equivalent initial state. """
    T = states.T[0]
    P = states.P[0]
    Y = states.Y[0].tolist()
    t = states.t[-1]
    n = states.shape[0]

    library = Path(__file__).resolve().parent

    code = dedent(f"""\
        addpath("{library}");
        graphics_toolkit("qt");
        GrafSimulation({T}, {P}, {Y}, {t}, {n}, "{saveas}");
        """)

    return code


def run_octave_program(states, saveas="graf_plot_octave"):
    """ Manage the code generation and execution of Octave code. """
    program_code = make_octave_program(states, saveas=saveas)

    with NamedTemporaryFile() as program_file:
        with open(program_file.name, "w") as fp:
            fp.write(program_code)

        run(["octave", program_file.name])

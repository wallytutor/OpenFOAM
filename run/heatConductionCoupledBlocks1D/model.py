# -*- coding: utf-8 -*-
import numpy as np
import pyvista as pv
from majordome_simulation.meshing import GmshOCCModel
from majordome_simulation.meshing import GeometricProgression
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
from scipy.special import erf
from screeninfo import get_monitors

RENDER: bool = False
""" Whether to render the Gmsh GUI during mesh generation. """

DEFAULT_OPTIONS: dict = {
    "General.Verbosity": 1,
    "Geometry.Points": True,
    "Geometry.Lines": True,
    "Geometry.Surfaces": True,
    "Mesh.ColorCarousel": 2,
    "Mesh.SurfaceFaces": True,
    "Mesh.SaveAll": False,
    "Mesh.SaveGroupsOfNodes": True,
    "Mesh.Algorithm": 6,
    "Mesh.ElementOrder": 1,
    "Mesh.MshFileVersion": 2.2
}
""" Default options for the Gmsh configuration. """

BLOCK: float = 4.0
""" Thickness of the first block (inches). """

DEPTH: float = 1.0
""" Total depth of the domain (large enough to approximate an infinite). """

GRID_FINE: float = 0.001
""" Grid size in first body and near the surface. """

GRID_COARSE: float = 0.1
""" Grid size at the bottom of the domain. """

RHO: float = 2500.0
""" Density of the material [kg/m³]. """

CP: float = 1000.0
""" Specific heat capacity of the material [J/(kg·K)]. """

K: float = 1.0
""" Thermal conductivity of the material [W/(m·K)]. """


def generate_domain(
        layer: float,
        saveas: str | None = None,
        overwrite: bool = False,
        options: dict = DEFAULT_OPTIONS
    ) -> None:
    """ Generates the mesh for the heat conduction 1D problem. """
    try:
        # Can this fail in HPC mode? "try block just in case"...
        # Monitor information for rendering the Gmsh GUI.
        monitor = get_monitors()[0]
        options["General.GraphicsWidth"]  = int(monitor.width * 0.8)
        options["General.GraphicsHeight"] = int(monitor.height * 0.8)
    except Exception:
        pass

    # Arbitrary, just for respecting OpenFOAM 3D requirements:
    side  = 0.1

    thickness1 = 0.0254 * layer
    thickness2 = DEPTH - thickness1

    n, q = GeometricProgression.fit(thickness2, GRID_FINE, GRID_COARSE)
    gp = GeometricProgression(n, GRID_FINE, q=q)

    layers1 = int(thickness1 / GRID_FINE)

    origin = (0.0, 0.0, 0.0)

    with GmshOCCModel(render=RENDER, **options) as model:
        #region: base
        tag_base = model.add_rectangle(*origin, side, side)
        model.synchronize()

        bounds = model.get_boundary([(2, tag_base)])
        bounds = [tag for dim, tag in bounds if dim == 1]

        for tag in bounds:
            model.set_transfinite_curve(tag, 2)

        model.set_recombine(2, tag_base)
        model.synchronize()
        #endregion: base

        #region: extrude1
        [tag_inter, vol1, *sides1] = model.extrude(
            dimTags     = [(2, tag_base)],
            dx          = 0.0,
            dy          = 0.0,
            dz          = thickness1,
            numElements = [layers1],
            recombine   = True
        )
        model.synchronize()

        tag_inter = tag_inter[1]
        vol1      = vol1[1]
        sides1    = [side[1] for side in sides1]
        #endregion: extrude1

        #region: extrude2
        # XXX this freezes, so I am brute-forcing the extrusion:
        # [tag_infty, vol2, *sides2] = model.extrude(
        #     dimTags     = [(2, tag_inter)],
        #     dx          = 0.0,
        #     dy          = 0.0,
        #     dz          = thickness2,
        #     numElements = [layers2],
        #     recombine   = True
        # )
        # model.synchronize()

        previous   = vol1
        to_extrude = (2, tag_inter)
        vol2   = []
        sides2 = []

        for dh in gp.heights:
            [tag_inter, vol, *sides] = model.extrude(
                dimTags     = [to_extrude],
                dx          = 0.0,
                dy          = 0.0,
                dz          = dh,
                numElements = [1],
                recombine   = True
            )
            model.synchronize()

            model.fragment([(3, previous)], [vol])
            model.synchronize()
            previous = vol[1]

            to_extrude = tag_inter
            vol2.append(vol[1])
            sides2.extend([side[1] for side in sides])
        #endregion: extrude2

        #region: physical groups

        surface     = [tag_base]
        sides_block = sides1
        sides_media = sides2
        infinite    = [to_extrude[1]]

        block = [vol1]
        media = vol2
        model.add_physical_groups(
            surfaces = [
                {"tags": surface,     "name": "surface",     "tag_id": 1},
                {"tags": sides_block, "name": "sides_block", "tag_id": 2},
                {"tags": sides_media, "name": "sides_media", "tag_id": 3},
                {"tags": infinite,    "name": "infinite",    "tag_id": 4},
            ],
            volumes = [
                {"tags": block, "name": "block", "tag_id": 10},
                {"tags": media, "name": "media", "tag_id": 20},
            ],
        )
        model.synchronize()
        #endregion: physical groups

        model.generate_mesh(dim=2)
        model.generate_mesh(dim=3)
        model.synchronize()

        if saveas and (overwrite or not Path(saveas).exists()):
            model.dump(saveas)


def kelvin_to_fahrenheit(T):
    """ Converts temperature from Kelvin to Fahrenheit. """
    return (T - 273.15) * 9/5 + 32


def analytical(z, t, Ti=297.04, Ts=810.97, alpha=K/(RHO*CP)):
    """ Semi-infinite medium analytical solution profile. """
    T = Ts + (Ti - Ts) * erf(-z / np.sqrt(4 * alpha * t))
    return kelvin_to_fahrenheit(T)


def get_profile(grid):
    """ Get temperature profile from an OpenFOAM solution."""
    def _get_profile(region):
        z = np.array(region.cell_centers().points[:, 2].tolist())
        T = np.array(region.cell_data["T"].tolist())
        return z, T

    block = grid["block"]["internalMesh"]
    media = grid["media"]["internalMesh"]

    z_block, T_block = _get_profile(block)
    z_media, T_media = _get_profile(media)

    z = np.concatenate((z_block, z_media))
    T = np.concatenate((T_block, T_media))
    return z, kelvin_to_fahrenheit(T)


def tabulate_profile(reader):
    """ Create a table with all time-step solutions. """
    matrix = []

    for this_time in reader.time_values[1:]:
        reader.set_active_time_value(this_time)
        matrix.append(get_profile(reader.read())[1])

    matrix = np.array(matrix)
    return matrix


def plot_evolution(z, times, data, z_max=0.05, saveas="evolution.png"):
    """ Plot the evolution of the temperature profile over time. """
    mask = z < z_max
    z_selected = z[mask]

    order = np.argsort(z_selected)
    z_selected = z_selected[order]

    selected = data[:, mask][:, order]
    times_s = np.array(times)
    times = times_s / 60

    z_lims = (z_selected[0], z_selected[-1])
    t_lims = (times[0], times[-1])
    extent = (z_lims[0], z_lims[1], t_lims[0], t_lims[1])

    Z, TT = np.meshgrid(z_selected, times_s)
    analytical_field = analytical(-Z, TT)

    plt.close("all")
    cb = plt.imshow(selected, aspect="auto", extent=extent,
                    origin="lower", cmap="hot")

    opts1 = dict(levels=[200.0], colors="cyan", linewidths=1.5)
    opts2 = dict(levels=[200.0], colors="lime", linewidths=1.5)

    line1 = Line2D([0], [0], color="cyan", linewidth=1.5,
                   label="Numerical 200 °F")
    line2 = Line2D([0], [0], color="lime", linewidth=1.5,
                   label="Analytical 200 °F")

    plt.contour(z_selected, times, selected, **opts1)
    plt.contour(z_selected, times, analytical_field, **opts2)
    plt.legend(handles=[line1, line2], loc=4, fontsize="small")

    plt.xlabel("z [m]")
    plt.ylabel("time [min]")
    plt.title("Temperature [°F]")

    cbar = plt.colorbar(cb, label="Temperature [°F]")
    cbar.set_ticks([100, 200, 400, 600, 800, 1000])
    plt.savefig(saveas, dpi=300)


def post_process():
    """ Post-process the case. """
    reader = pv.POpenFOAMReader("case.foam")
    reader.set_active_time_value(reader.time_values[1])

    z, _  = get_profile(reader.read())
    data  = tabulate_profile(reader)
    times = reader.time_values[1:]

    plot_evolution(z, times, data)


if __name__ == "__main__":
    if not Path("mesh.msh").exists():
        generate_domain(BLOCK, saveas="mesh.msh", overwrite=False)
    else:
        post_process()

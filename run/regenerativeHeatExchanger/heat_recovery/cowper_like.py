# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pyvista as pv
import majordome_simulation.meshing as ms
from foamlib import FoamFile
from pathlib import Path
from majordome_engineering.transport import SolutionDimless
from majordome_simulation.meshing import GmshOCCModel
from majordome_simulation.meshing import RingBuilder
from majordome_simulation.meshing import CircularCrossSection
from majordome_utilities.plotting import plot2d
from ruamel.yaml import YAML
from screeninfo import get_monitors

from .calculators import SkinFrictionFactor
from .calculators import WallGradingCalculator
from .thermocline import ThermoclineModel

MONITOR = get_monitors()[0]

DEFAULT_OPTIONS = {
    "General.GraphicsWidth": int(MONITOR.width*0.8),
    "General.GraphicsHeight": int(MONITOR.height*0.8),
    "General.Verbosity": 5,
    "Geometry.Points": True,
    "Geometry.Lines": True,
    "Geometry.Surfaces": True,
    "Mesh.CharacteristicLengthMin": 1.0e-05,
    "Mesh.CharacteristicLengthMax": 1.0e-02,
    "Mesh.ColorCarousel": 2,
    "Mesh.SurfaceFaces": True,
    "Mesh.SaveAll": False,
    "Mesh.SaveGroupsOfNodes": True,
    "Mesh.Algorithm": 6,
    "Mesh.ElementOrder": 1,
    "Mesh.MshFileVersion": 2.2
}


class CowperLikeGeometry:
    """ Geometry of a Cowper-like regenerative heat exchanger.

    It is represented by a fluid hole surrounded by a hexagonal solid
    region, extruded in the z-direction. The geometry is created using
    the Gmsh Python API and meshed with transfinite meshing in the radial
    direction and extrusion in the z-direction.

    Parameters
    ----------
    m_h: float
        Ratio of the distance from the center of the fluid hole to the
        face of the hexagon, to the diameter of the fluid hole.
    D_h: float
        Diameter of the fluid hole.
    h_t: float
        Height of the extruded geometry in the z-direction.
    num_points_angular: int
        Number of points along the angular direction of the fluid hole
        and hexagon for mesh structuring.
    fluid_bl_tot: float
        Total thickness of the boundary layer in the fluid region.
    fluid_bl_ext: float
        Thickness of the outermost boundary layer cell in the fluid region.
    fluid_bl_int: float
        Thickness of the innermost boundary layer cell in the fluid region.
    solid_bl_ext: float
        Thickness of the outermost boundary layer cell in the solid region.
    solid_bl_int: float
        Thickness of the innermost boundary layer cell in the solid region.
    rel_layer: float
        Relative thickness of each layer in the extrusion direction, as a
        fraction of the diameter of the fluid hole.
    num_points_core: int
        Number of points along the core line in the radial direction for
        mesh structuring.
    """
    def __init__(self, *,
            m_h: float,
            D_h: float,
            h_t: float,
            num_points_angular: int = 6,
            fluid_bl_tot: float = 0.003,
            fluid_bl_ext: float = 0.0001,
            fluid_bl_int: float = 0.0006,
            solid_bl_ext: float = 0.0050,
            solid_bl_int: float = 0.0001,
            rel_layer: float    = 0.15,
            core_radius_fraction: float = 0.75,
            num_points_core: int = 2,
        ) -> None:
        self.m_h = m_h
        self.D_h = D_h
        self.h_t = h_t
        self.num_points_angular = num_points_angular
        self.fluid_bl_tot = fluid_bl_tot
        self.fluid_bl_ext = fluid_bl_ext
        self.fluid_bl_int = fluid_bl_int
        self.solid_bl_ext = solid_bl_ext
        self.solid_bl_int = solid_bl_int
        self.core_radius_fraction = core_radius_fraction
        self.num_points_core = num_points_core

        # Six-fold symmetry of the system:
        self.num_splits = 6

        # Rotate so that face of the hexagon is aligned with the x-axis:
        self.rotation   = 30

        # Number of layers in the extrusion direction:
        self.n_layers   = int(h_t / (rel_layer * D_h))

        # Apothem of the hexagon, i.e. distance from center to face:
        self.center_to_wall = (1 + m_h) * D_h / 2

        # Outer radius of the hexagon and fluid hole:
        self.R_out_solid = 2 * self.center_to_wall / np.sqrt(3)
        self.R_out_fluid = D_h / 2

    def _create_fluid(self, model):
        fluid_hole = CircularCrossSection(
            model              = model,
            radius             = self.R_out_fluid,
            boundary_thickness = self.fluid_bl_tot,
            cell_size_external = self.fluid_bl_ext,
            cell_size_internal = self.fluid_bl_int,
            num_points_angular = self.num_points_angular,
            num_splits         = self.num_splits,
            core_unstructured  = False,
            core_polygonal     = True,
            radius_fraction    = self.core_radius_fraction,
            rotation           = self.rotation
        )
        return fluid_hole

    def _create_hexagon(self, model, fluid_hole):
        solid_bl_tot = self.R_out_solid - self.R_out_fluid

        callbacks = RingBuilder.get_progression_callbacks(
            model,
            self.num_points_angular,
            solid_bl_tot,
            self.solid_bl_ext,
            self.solid_bl_int
        )

        callback_lines, callback_surfaces = callbacks

        hexa = RingBuilder(
            model              = model,
            splits             = self.num_splits,
            radius_out         = self.R_out_solid,
            points_in          = fluid_hole.ring.points_out,
            lines_in           = fluid_hole.ring.lines_out,
            callback_lines     = callback_lines,
            callback_surfaces  = callback_surfaces,
            rotation           = self.rotation
        )
        return hexa

    def _create_volume(self, model, fluid_hole, hexa):
        hole, core = fluid_hole.ring, fluid_hole.core

        # TODO add API to access this:
        p_orig = fluid_hole._p_origin

        p0 = core.points_in[-1]
        p1 = core.points_in[0]

        l0 = core.lines_in[-1]
        l1 = model.add_line(p1, p_orig)
        l2 = model.add_line(p_orig, p0)

        loop = model.add_curve_loop([l0, l1, l2])
        surf = model.add_plane_surface([loop])
        model.synchronize()

        base = [(2, surf),
                (2, core.surface_tags[-1]),
                (2, hole.surface_tags[-1]),
                (2, hexa.surface_tags[-1])]

        removes  = core.surface_tags[:-1]
        removes += hole.surface_tags[:-1]
        removes += hexa.surface_tags[:-1]
        model.remove([(2, item) for item in removes])
        model.synchronize()

        model.set_transfinite_curve(l1, self.num_points_core)
        model.set_transfinite_curve(l2, self.num_points_core)

        opts = dict(numElements=[self.n_layers], recombine=True)
        tags = model.extrude(base, 0, 0, self.h_t, **opts)
        return tags, base

    def _create_labels(self, model, tags, base):
        extruded = ms.get_extrusion_tags(tags, 2)
        extruded_super, extruded_ndim = extruded

        fluid      = extruded_super[:3]
        outlet     = extruded_ndim[:3]
        solid      = [extruded_super[3]]
        top        = [extruded_ndim[3]]
        inlet      = [s[1] for s in base[:3]]
        bottom     = [base[3][1]]
        symmetry   = [34]
        fluidLeft  = [21, 25, 29]
        fluidRight = [22, 24, 28]
        solidLeft  = [33]
        solidRight = [32]

        model.add_physical_groups(
            surfaces=[
                {"tags": inlet,      "name": "fluidInlet",    "tag_id": 10},
                {"tags": outlet,     "name": "fluidOutlet",   "tag_id": 20},
                {"tags": bottom,     "name": "solidBottom",   "tag_id": 30},
                {"tags": top,        "name": "solidTop",      "tag_id": 40},
                {"tags": symmetry,   "name": "solidSymmetry", "tag_id": 50},
                {"tags": fluidLeft,  "name": "fluidLeft",     "tag_id": 60},
                {"tags": fluidRight, "name": "fluidRight",    "tag_id": 61},
                {"tags": solidLeft,  "name": "solidLeft",     "tag_id": 70},
                {"tags": solidRight, "name": "solidRight",    "tag_id": 71},
            ],
            volumes=[
                {"tags": fluid, "name": "fluid", "tag_id": 100},
                {"tags": solid, "name": "solid", "tag_id": 200},
            ]
        )

    def create_model(self,
            saveas: str | None = None,
            render: bool=True,
            options: dict = DEFAULT_OPTIONS
        ) -> None:
        with GmshOCCModel(render=render, **options) as model:
            fluid_hole = self._create_fluid(model)
            model.synchronize()

            hexa = self._create_hexagon(model, fluid_hole)
            model.synchronize()

            tags, base = self._create_volume(model, fluid_hole, hexa)
            model.synchronize()

            self._create_labels(model, tags, base)
            model.synchronize()

            model.generate_mesh(dim=2)
            model.generate_mesh(dim=3)
            model.synchronize()

            if saveas:
                model.dump(saveas)


def get_model_dump_mesh(config="dimensioning.yaml", mesh="mesh.msh"):
    model = ThermoclineModel(config)
    sol = SolutionDimless("airish.yaml")

    y_p = model.num_y_plus
    U_h = model.fn_U_g(*model.args_charging)

    sol.set_state(model.num_T_h, 101325, "N2: 0.79, O2: 0.21")
    calc = WallGradingCalculator.from_solution(sol, L=model.num_D_h, U=U_h)

    def f_tur(Re):
        return SkinFrictionFactor.smooth_wall(Re, check=False) / 8

    y_first = calc.first_layer(y_p, skin_factor=f_tur)

    yaml = YAML()
    data = yaml.load(open("dimensioning.yaml"))

    if not Path(mesh).exists():
        geom = CowperLikeGeometry(
            m_h = data["m_h"],
            D_h = data["D_h"],
            h_t = data["h_t"],

            # Use with care based on the above values:
            num_points_angular = 6,
            num_points_core    = 4,
            core_radius_fraction = 0.8,
            fluid_bl_tot       = 0.005,
            fluid_bl_ext       = 0.5*y_first,
            fluid_bl_int       = 1.5*y_first,
            solid_bl_ext       = 3.0*y_first,
            solid_bl_int       = 0.5*y_first,
            rel_layer          = 0.25,
        )

        # mesh = None # DEBUG (no write)
        geom.create_model(saveas=mesh, render=True)

    return model


def make_charging(model):
    """ Create the charging inputs file.

    Notice:
    - The inlet velocity is negative because the fluid enters from the
      *outlet* side (top) of the geometry, poiting downwards.
    """
    with FoamFile("constant/userParameters") as f:
        f["fluidTemperature"] = model.num_T_h
        f["solidTemperature"] = model.num_T_c
        f["inletVelocity"] = -1.0 * model.fn_U_g(*model.args_charging)
        f["duration"] = model.num_H_t / model.num_P_c
        f["nProcs"] = r"$NUM_PROCS"

    def add_property(f, name, value):
        f[name] = {}
        f[name]["type"] = "uniform"
        f[name]["value"] = value

    with FoamFile("constant/solid/physicalProperties") as f:
        f["thermoType"] = "constSolidThermo"
        add_property(f, "rho",   model.num_rho_s)
        add_property(f, "Cv",    model.num_c_ps)
        add_property(f, "kappa", model.num_k_s)


#region: tabular post-processing
def get_post(region, function, fname, time="0.000000e+00"):
    """ Get the path to a post-processing table. """
    return "/".join(["postProcessing", region, function, time, f"{fname}.dat"])


def load_table(region, function, fname, time="0.000000e+00"):
    """ Load a table from the post-processing directory. """
    fname = get_post(region, function, fname, time)
    return pd.read_csv(fname, sep=r"\s+", header=None, comment="#")


def plot_table(region, function, fname, time="0.000000e+00", ylabel=None):
    """ Quick plotting of a table from the post-processing directory.

    This is for debugging, please consider using the next functions.
    """
    df = load_table(region, function, fname, time)

    if function.startswith("pressure"):
        df[1] -= 101325.0

    p = plot2d(df[0], df[1])
    p.axes[0].set_xlabel("Time [s]")
    p.axes[0].set_ylabel(ylabel or function)
    return p


def plot_temperature(title, origin="0.000000e+00", **kws):
    """ Plot the temperature from the post-processing directory. """
    def loader(fname, kind="surfaceFieldValue"):
        return load_table("fluid", fname,  kind, time=origin)

    dfi = loader("temperatureInlet")
    dfo = loader("temperatureOutlet")

    t1 = dfi[0].to_numpy()
    t2 = dfo[0].to_numpy()
    p1 = dfi[1].to_numpy() - 273.15
    p2 = dfo[1].to_numpy() - 273.15

    p = plot2d(t1, p1, color="b", label="Cold side")
    p.axes[0].plot(t2, p2, color="r", label="Hot side")
    p.axes[0].set_title(title)
    p.axes[0].set_xlabel("Time [s]")
    p.axes[0].set_ylabel("Temperature [°C]")
    p.axes[0].legend(loc=kws.get("loc", 4), fontsize="small")
    return p


def plot_pressure(title, origin="0.000000e+00", **kws):
    """ Plot the pressure from the post-processing directory. """
    def loader(fname, kind="surfaceFieldValue"):
        return load_table("fluid", fname,  kind, time=origin)

    dfi = loader("pressureInlet")
    dfo = loader("pressureOutlet")

    t1 = dfi[0].to_numpy()
    t2 = dfo[0].to_numpy()
    p1 = (dfi[1].to_numpy() - 101325.0) / 100.0
    p2 = (dfo[1].to_numpy() - 101325.0) / 100.0

    p = plot2d(t1, p1, color="b", label="Cold side")
    p.axes[0].plot(t2, p2, color="r", label="Hot side")
    p.axes[0].set_title(title)
    p.axes[0].set_xlabel("Time [s]")
    p.axes[0].set_ylabel("Pressure [mbar]")
    p.axes[0].legend(loc=kws.get("loc", 4), fontsize="small")
    return p


def plot_flowrate(title, origin="0.000000e+00", **kws):
    """ Plot the flow rate from the post-processing directory. """
    def loader(fname, kind="surfaceFieldValue"):
        return load_table("fluid", fname,  kind, time=origin)

    dfi = loader("flowRateInlet")
    dfo = loader("flowRateOutlet")

    t1 = dfi[0].to_numpy()
    t2 = dfo[0].to_numpy()
    p1 = dfi[1].to_numpy() * 6 * 1000
    p2 = dfo[1].to_numpy() * 6 * 1000

    s = min(len(t1), len(t2))
    tb = t1[:s]
    pb = p1[:s] + p2[:s]

    p = plot2d(t1, p1, color="b", label="Cold side")
    p.axes[0].plot(t2, p2, color="r", label="Hot side")
    p.axes[0].plot(tb, pb, color="k", label="Balance")
    p.axes[0].set_title(title)
    p.axes[0].set_xlabel("Time [s]")
    p.axes[0].set_ylabel("Flow rate [g/s]")
    p.axes[0].legend(loc=kws.get("loc", 4), fontsize="small")
    return p


def plot_convergence():
    """ Plot the convergence from the log files. """
    final = ["UxFinalRes_0", "UyFinalRes_0", "UzFinalRes_0",
             "eFinalRes_0", "hFinalRes_0",
             "kFinalRes_0", "omegaFinalRes_0"]

    p = plot2d(0, 0)

    for fname in final:
        try:
            df = pd.read_csv(f"logs/{fname}", sep=r"\s+", header=None)
            p.axes[0].plot(np.log10(df[1]), label=fname.split("_")[0])
        except FileNotFoundError as err:
            print(err)

    p.axes[0].set_xlabel("Time step")
    p.axes[0].set_ylabel("log10(residual)")
    p.axes[0].legend(loc=1, fontsize="small", ncol=2)
#endregion: tabular post-processing

#region: graphical post-processing
class CowperLikePost:
    __slots__ = (
        "_reader",
        "_scale",
        "_time",
        "_fluid_internal",
        "_solid_internal",
        "_slice_fluid",
        "_slice_solid",
    )

    def __init__(self, scale=(1, 1, 1)):
        self._reader = pv.POpenFOAMReader("case.foam")
        self._scale = scale

    def _get_slice(self, internal_mesh):
        data = internal_mesh.slice("y")

        if "U" in data.cell_data:
            Umag = np.linalg.norm(data.cell_data["U"], axis=1)
            data.cell_data["Umag"] = Umag

        return data.scale(self._scale, inplace=False)

    @property
    def time_values(self):
        return self._reader.time_values

    @property
    def slice_fluid(self):
        return self._slice_fluid.copy()

    @property
    def slice_solid(self):
        return self._slice_solid.copy()

    @staticmethod
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

    @staticmethod
    def align_camera(plot, xc=0.025, zc=0.02, ps=0.017):
        plot.renderer.set_background("#FFFFFF")
        plot.camera_position = CowperLikePost.camera_position_xz(xc, zc)
        plot.camera.parallel_projection = True
        plot.camera.parallel_scale = ps

    def load_state(self, time=None):
        if time is None:
            time = self.time_values[-1]

        self._time = time
        self._reader.set_active_time_value(time)

        mesh = self._reader.read()
        self._fluid_internal = mesh["fluid"]["internalMesh"]
        self._solid_internal = mesh["solid"]["internalMesh"]

        self._slice_fluid = self._get_slice(self._fluid_internal)
        self._slice_solid = self._get_slice(self._solid_internal)

        pRel = self._slice_fluid.cell_data["p"] - 101325.0
        self._slice_fluid.cell_data["pRel"] = pRel

    def get_options(self, **kws) -> dict:
        fn_title = kws.pop("fn_title", lambda t: f"At t = {t:.0f} s\n")

        opts = {
            "pbr": False,
            "show_edges": False,
            "scalar_bar_args": {
                "title": fn_title(self._time),
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
        }
        opts.update(kws)
        return opts
#endregion: graphical post-processing
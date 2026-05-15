# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pyvista as pv
import warnings
import majordome.simulation as ms
from foamlib import FoamFile
from majordome import constants
from majordome.simulation import FoamPostProcessingLoader
from majordome.utilities import plot_xy
from screeninfo import get_monitors

from .thermocline import ThermoclineModel

DEFAULT_OPTIONS = {
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


def _simple_warning(message, category, filename, lineno, file=None, line=None):
    return f"{category.__name__}: {message}\n"


warnings.formatwarning = _simple_warning


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
        fluid_hole = ms.CircularCrossSection(
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

        callbacks = ms.RingBuilder.get_progression_callbacks(
            model,
            self.num_points_angular,
            solid_bl_tot,
            self.solid_bl_ext,
            self.solid_bl_int
        )

        callback_lines, callback_surfaces = callbacks

        hexa = ms.RingBuilder(
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
                {"tags": inlet,      "name": "fluidCold",    "tag_id": 10},
                {"tags": outlet,     "name": "fluidHot",   "tag_id": 20},
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
            render: bool = True,
            options: dict = DEFAULT_OPTIONS
        ) -> None:

        try:
            # Can this fail in HPC mode? "try block just in case"...
            # Monitor information for rendering the Gmsh GUI.
            monitor = get_monitors()[0]
            options["General.GraphicsWidth"]  = int(monitor.width * 0.8)
            options["General.GraphicsHeight"] = int(monitor.height * 0.8)
        except Exception:
            pass

        with ms.GmshOCCModel(render=render, **options) as model:
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


def make_solid(model):
    """ Create the solid physical properties file. """
    def add_property(f, name, value):
        f[name] = {}
        f[name]["type"] = "uniform"
        f[name]["value"] = value

    with FoamFile("constant/solid/physicalProperties") as f:
        f["thermoType"] = "constSolidThermo"
        add_property(f, "rho",   model.num_rho_s)
        add_property(f, "Cv",    model.num_c_ps)
        add_property(f, "kappa", model.num_k_s)


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


class PostProcessing:
    """ Unified post-processing for the heat exchanger case. """
    def __init__(self, root, mode, domain=None, **kw):
        if domain is None:
            self.fluid = FoamPostProcessingLoader(domain=None, root=root)
            self.solid = None
        else:
            self.fluid = FoamPostProcessingLoader(domain="fluid", root=root)
            self.solid = FoamPostProcessingLoader(domain="solid", root=root)

        self.case = root
        self.root = self.fluid.root_directory
        self.mode = mode.lower()
        self.logs = self.root.parent / "logs"

        self.auto_load_fields = kw.get("auto_load_fields", True)
        self.time   = None
        self.fields = {}

        self._load_fields(**kw)

    #region: field loading
    def _load_fields(self, **kw):
        """ Load the fields for post-processing. """
        try:
            PostProcessing._ensure_foam_case(self.case)

            model = ThermoclineModel(kw.get("config", "dimensioning.yaml"))
            scale = (1, 1, model.num_D_h / model.num_h_t)

            reader = pv.POpenFOAMReader(self.case / "case.foam")
            self.time = reader.time_values[-1]
            reader.set_active_time_value(self.time)

            mesh = reader.read()

            if self.solid is None:
                fluid = PostProcessing._load_region(mesh, scale, name=None)
                self.fields = {"fluid": fluid}
            else:
                fluid = PostProcessing._load_region(mesh, scale, name="fluid")
                solid = PostProcessing._load_region(mesh, scale, name="solid")
                self.fields = {"fluid": fluid, "solid": solid}
        except Exception as err:
            warnings.warn(f"Failed to load fields for post-processing: {err}")

    @staticmethod
    def _ensure_foam_case(case_dir):
        if not (case_dir / "case.foam").exists():
            with open(case_dir / "case.foam", "w") as f:
                f.write("dummy")

    @staticmethod
    def _load_region(mesh, scale, name=None):
        if name is None:
            data = mesh.slice("y")[0]
            data.points *= scale
        else:
            data = mesh[name].slice("y")[0]
            data.points *= scale

        if "rho" in data.cell_data:
            data.cell_data["rho"] = data.cell_data["rho"] * 1000.0

        if "p" in data.cell_data:
            data.cell_data["p"] = data.cell_data["p"] - constants.P_NORMAL

        if "T" in data.cell_data:
            data.cell_data["T"] = data.cell_data["T"] - constants.T_NORMAL

        return data.cell_data_to_point_data()

    @staticmethod
    def _ensure_has_time(fn):
        """ Decorator to ensure that time is available for plotting. """
        def wrapper(self, *args, **kwargs):
            if self.auto_load_fields:
                self._load_fields(**kwargs)

            if self.time is None or not self.fields:
                warnings.warn("Time information is not available for plotting.")
                return
            return fn(self, *args, **kwargs)
        return wrapper
    #endregion: field loading

    #region: postProcessing
    def plot_inlet_mass_flow_rate(self):
        """ Plot the inlet mass flow rate over time. """
        df1 = self.fluid.load_report("flowRateHot")
        df2 = self.fluid.load_report("flowRateCold")

        n = min(len(df1), len(df2))

        x = df1.iloc[:n, 0]
        m1 = df1.iloc[:n, 1] * (-1000)
        m2 = df2.iloc[:n, 1] * (-1000)

        p = plot_xy()
        p.axes[0].plot(x, m1, label="Hot Side")
        p.axes[0].plot(x, m2, label="Cold Side")
        p.axes[0].legend(loc="best", fontsize="small")
        p.axes[0].set_title(f"Inlet mass flow rate during {self.mode}")
        p.axes[0].set_xlabel("Time [s]")
        p.axes[0].set_ylabel("Mass flow rate [g/s]")
        p.resize(6, 5)
        p.savefig(self.root / "01-inlet_mass_flow_rate.png")

    def plot_imbalance_mass_flow_rate(self):
        """ Plot the imbalance of mass flow rate over time. """
        df1 = self.fluid.load_report("flowRateHot")
        df2 = self.fluid.load_report("flowRateCold")

        n = min(len(df1), len(df2))

        x = df1.iloc[:n, 0]
        m1 = df1.iloc[:n, 1] * (-1000)
        m2 = df2.iloc[:n, 1] * (-1000)
        imbalance = (m1 + m2) / m1 * 100

        p = plot_xy(x, imbalance)
        p.axes[0].set_title(f"Mass flow imbalance during {self.mode}")
        p.axes[0].set_xlabel("Time [s]")
        p.axes[0].set_ylabel("Mass flow imbalance [%]")
        p.resize(6, 5)
        p.savefig(self.root / "02-imbalance_mass_flow_rate.png")

    def plot_total_pressure(self):
        """ Plot the total pressure over time. """
        df1 = self.fluid.load_report("pressureCold")
        df2 = self.fluid.load_report("pressureHot")

        n = min(len(df1), len(df2))

        x = df1.iloc[:n, 0]
        p1 = (df1.iloc[:n, 1] - constants.P_NORMAL) / 100
        p2 = (df2.iloc[:n, 1] - constants.P_NORMAL) / 100

        p = plot_xy()
        p.axes[0].plot(x, p1, label="Cold Side")
        p.axes[0].plot(x, p2, label="Hot Side")
        p.axes[0].legend(loc="best", fontsize="small")
        p.axes[0].set_title(f"Total pressure during {self.mode}")
        p.axes[0].set_xlabel("Time [s]")
        p.axes[0].set_ylabel("Total pressure [mbar]")
        p.resize(6, 5)
        p.savefig(self.root / "03-total_pressure.png")

    def plot_pressure_drop(self):
        """ Plot the pressure drop over time. """
        df1 = self.fluid.load_report("pressureCold")
        df2 = self.fluid.load_report("pressureHot")

        n = min(len(df1), len(df2))

        x = df1.iloc[:n, 0]
        p1 = (df1.iloc[:n, 1] - constants.P_NORMAL) / 100
        p2 = (df2.iloc[:n, 1] - constants.P_NORMAL) / 100
        dp = p2 - p1 if self.mode == "charging" else p1 - p2

        p = plot_xy(x, dp)
        p.axes[0].set_title(f"Pressure drop during {self.mode}")
        p.axes[0].set_xlabel("Time [s]")
        p.axes[0].set_ylabel("Pressure drop [mbar]")
        p.resize(6, 5)
        p.savefig(self.root / "04-pressure_drop.png")

    def plot_total_temperature(self):
        """ Plot the total temperature over time. """
        df1 = self.fluid.load_report("temperatureCold")
        df2 = self.fluid.load_report("temperatureHot")

        n = min(len(df1), len(df2))

        x = df1.iloc[:n, 0]
        p1 = df1.iloc[:n, 1] - constants.T_NORMAL
        p2 = df2.iloc[:n, 1] - constants.T_NORMAL

        p = plot_xy()
        p.axes[0].plot(x, p1, label="Cold Side")
        p.axes[0].plot(x, p2, label="Hot Side")
        p.axes[0].legend(loc="best", fontsize="small")
        p.axes[0].set_title(f"Total temperature during {self.mode}")
        p.axes[0].set_xlabel("Time [s]")
        p.axes[0].set_ylabel("Total temperature [°C]")
        p.resize(6, 5)
        p.savefig(self.root / "05-total_temperature.png")

    def plot_y_plus(self):
        """ Plot the y+ values over time. """
        df = self.fluid.load_report("yPlus")

        x = df.iloc[:, 0]
        p_min = df.iloc[:, 2]
        p_max = df.iloc[:, 3]
        p_avg = df.iloc[:, 4]

        p = plot_xy()
        p.axes[0].plot(x, p_min, label="Min y+")
        p.axes[0].plot(x, p_max, label="Max y+")
        p.axes[0].plot(x, p_avg, label="Avg y+")
        p.axes[0].legend(loc="best", fontsize="small")
        p.axes[0].set_title(f"y+ values during {self.mode}")
        p.axes[0].set_xlabel("Time [s]")
        p.axes[0].set_ylabel("y+")
        p.resize(6, 5)
        p.savefig(self.root / "06-y_plus.png")
    #endregion: postProcessing

    #region: foamLog
    def plot_convergence(self, ymin=-10, ymax=-6):
        """ Plot the convergence from the log files. """
        final = ["UxFinalRes_0", "UyFinalRes_0", "UzFinalRes_0",
                 "eFinalRes_0", "hFinalRes_0",
                 "kFinalRes_0", "omegaFinalRes_0"]

        has_plots = False
        p = plot_xy()

        for fname in final:
            try:
                df = pd.read_csv(self.logs / fname, sep=r"\s+", header=None)
                p.axes[0].plot(np.log10(df[1]), label=fname.split("_")[0])
                has_plots = True
            except FileNotFoundError as err:
                print(err)

        if not has_plots:
            warnings.warn("No convergence data found to plot.")
            return

        p.axes[0].set_ylim(ymin, ymax)
        p.axes[0].set_title(f"Convergence during {self.mode}")
        p.axes[0].set_xlabel("Time step")
        p.axes[0].set_ylabel("log10(residual)")
        p.axes[0].legend(loc=1, fontsize="small", ncol=2)
        p.resize(6, 5)
        p.savefig(self.root / "convergence.png")
    #endregion: foamLog

    #region: fields
    @staticmethod
    def get_options(time, scalar_bar={}, **kws) -> dict:
        fn_title = kws.pop("fn_title", lambda t: f"At t = {t:.0f} s\n")

        scalar_bar_args = {
            "title": fn_title(time),
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
        scalar_bar_args.update(scalar_bar)

        opts = {
            "pbr": False,
            "show_edges": False,
            "scalar_bar_args": scalar_bar_args,
        }
        opts.update(kws)

        return opts

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
        plot.camera_position = PostProcessing.camera_position_xz(xc, zc)
        plot.camera.parallel_projection = True
        plot.camera.parallel_scale = ps

    @_ensure_has_time
    def plot_field_temperature(self, show_solid=False, show_fluid=True):
        if not show_solid and not show_fluid:
            warnings.warn("At least one of domains must be displayed.")
            return

        if (fluid := self.fields.get("fluid", None)) is None:
            warnings.warn("Fluid data not found.")
            return

        if "T" not in fluid.point_data:
            warnings.warn("Temperature field 'T' not found in the data.")
            return

        solid = None

        if show_solid and self.solid is not None and "solid" in self.fields:
            solid = self.fields["solid"]

        opts = self.get_options(
            self.time, scalar_bar={"fmt": "%.0f"},
            fn_title=lambda t: f"Temperature at t = {t:.1f} s\n"
        )

        plot = pv.Plotter()

        if show_fluid:
            plot.add_mesh(fluid, scalars="T", cmap="fire", **opts)

        if show_solid and solid is not None:
            plot.add_mesh(solid, scalars="T", cmap="fire", **opts)

        if show_fluid and not show_solid:
            self.align_camera(plot, xc=0.0125, zc=0.02, ps=0.017)
        elif show_solid and show_fluid:
            self.align_camera(plot, xc=0.0250, zc=0.02, ps=0.017)
        else:
            self.align_camera(plot, xc=0.0375, zc=0.02, ps=0.017)

        plot.show()

    @_ensure_has_time
    def plot_field_pressure(self):
        if (fluid := self.fields.get("fluid", None)) is None:
            warnings.warn("Fluid data not found.")
            return

        if "p" not in fluid.point_data:
            warnings.warn("Pressure field 'p' not found in the data.")
            return

        opts = self.get_options(
            self.time, scalar_bar={"fmt": "%.0f"},
            fn_title=lambda t: f"Pressure at t = {t:.1f} s\n"
        )

        plot = pv.Plotter()
        plot.add_mesh(fluid, scalars="p", cmap="plasma", **opts)
        self.align_camera(plot, xc=0.0125, zc=0.02, ps=0.017)
        plot.show()

    @_ensure_has_time
    def plot_field_density(self):
        if (fluid := self.fields.get("fluid", None)) is None:
            warnings.warn("Fluid data not found.")
            return

        if "rho" not in fluid.point_data:
            warnings.warn("Density field 'rho' not found in the data.")
            return

        opts = self.get_options(
            self.time, scalar_bar={"fmt": "%.2f"},
            fn_title=lambda t: f"Density at t = {t:.1f} s\n"
        )

        plot = pv.Plotter()
        plot.add_mesh(fluid, scalars="rho", cmap="plasma", **opts)
        self.align_camera(plot, xc=0.0125, zc=0.02, ps=0.017)
        plot.show()

    @_ensure_has_time
    def plot_field_velocity(self):
        if (fluid := self.fields.get("fluid", None)) is None:
            warnings.warn("Fluid data not found.")
            return

        if "U" not in fluid.point_data:
            warnings.warn("Velocity field 'U' not found in the data.")
            return

        opts = self.get_options(
            self.time, scalar_bar={"fmt": "%.2f"},
            fn_title=lambda t: f"Velocity at t = {t:.1f} s\n"
        )

        plot = pv.Plotter()
        plot.add_mesh(fluid, scalars="U", cmap="jet", **opts)
        self.align_camera(plot, xc=0.0125, zc=0.02, ps=0.017)
        plot.show()

    @_ensure_has_time
    def plot_field_buoyancy_pressure(self):
        if (fluid := self.fields.get("fluid", None)) is None:
            warnings.warn("Fluid data not found.")
            return

        if "p_rgh" not in fluid.point_data:
            warnings.warn("Pressure field 'p_rgh' not found in the data.")
            return

        opts = self.get_options(
            self.time, scalar_bar={"fmt": "%.0f"},
            fn_title=lambda t: f"Buoyancy pressure at t = {t:.1f} s\n"
        )

        plot = pv.Plotter()
        plot.add_mesh(fluid, scalars="p_rgh", cmap="plasma", **opts)
        self.align_camera(plot, xc=0.0125, zc=0.02, ps=0.017)
        plot.show()
    #endregion: fields

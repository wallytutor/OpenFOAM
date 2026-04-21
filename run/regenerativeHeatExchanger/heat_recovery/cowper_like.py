# -*- coding: utf-8 -*-
import numpy as np
import majordome_simulation.meshing as ms
from foamlib import FoamFile
from majordome_simulation.meshing import GmshOCCModel
from majordome_simulation.meshing import RingBuilder
from majordome_simulation.meshing import CircularCrossSection

from screeninfo import get_monitors

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
            radius_fraction    = 0.5,
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

        model.set_transfinite_curve(l1, 5)
        model.set_transfinite_curve(l2, 5)

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


def make_charging(model):
    with FoamFile("constant/userParameters") as f:
        f["fluidTemperature"] = model.num_T_h
        f["solidTemperature"] = model.num_T_c
        f["inletVelocity"] = model.fn_U_g(*model.args_charging)
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

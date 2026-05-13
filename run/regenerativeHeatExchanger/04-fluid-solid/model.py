# -*- coding: utf-8 -*-
import os
from pathlib import Path
from heat_recovery.cowper_like import CowperLikeGeometry
from heat_recovery.cowper_like import PostProcessing
from heat_recovery.cowper_like import make_solid
from heat_recovery.cowper_like import make_charging
from heat_recovery.thermocline import ThermoclineModel
from majordome.engineering import SolutionDimless
from majordome.engineering import SkinFrictionFactor
from majordome.engineering import WallGradingCalculator
from ruamel.yaml import YAML

RENDER: bool = False
""" Whether to render the Gmsh GUI during mesh generation. """


def generate_charging_case(config, mesh) -> None:
    """ Generate the mesh and conditions for the charging case. """
    model = ThermoclineModel(config)
    # make_solid(model)
    # make_charging(model)

    y_p = model.num_y_plus
    U_h = model.fn_U_g(*model.args_charging)

    sol = SolutionDimless("airish.yaml")
    sol.set_state(model.num_T_h, 101325, "N2: 0.79, O2: 0.21")
    calc = WallGradingCalculator.from_solution(sol, L=model.num_D_h, U=U_h)

    def f_tur(Re):
        return SkinFrictionFactor.smooth_wall(Re, check=False) / 8

    y_first = calc.first_layer(y_p, skin_factor=f_tur)

    if not Path(mesh).exists():
        data = YAML().load(open(config))

        geom = CowperLikeGeometry(
            m_h = data["m_h"],
            D_h = data["D_h"],
            h_t = data["h_t"],

            # Use with care based on the above values:
            num_points_angular = 6,
            num_points_core    = 4,
            core_radius_fraction = 0.8,
            fluid_bl_tot       = 0.005,
            fluid_bl_ext       = 0.5 * y_first,
            fluid_bl_int       = 1.5 * y_first,
            solid_bl_ext       = 3.0 * y_first,
            solid_bl_int       = 0.5 * y_first,
            rel_layer          = 0.25,
        )

        # mesh = None # DEBUG (no write)
        geom.create_model(saveas=mesh, render=RENDER)


def post_process(here):
    post = PostProcessing(here, mode="charging", domain=None)

    post.plot_inlet_mass_flow_rate()
    post.plot_imbalance_mass_flow_rate()
    post.plot_total_pressure()
    post.plot_pressure_drop()
    post.plot_total_temperature()
    post.plot_y_plus()
    post.plot_convergence()

    post.plot_field_temperature()
    post.plot_field_pressure()
    post.plot_field_density()
    post.plot_field_velocity()
    # post.plot_field_turb_kinetic_energy()
    # post.plot_field_turb_dissipation_rate()
    post.plot_field_buoyancy_pressure()


if __name__ == "__main__":
    try:
        here = Path(__file__).parent
    except NameError:
        here = Path(os.getcwd()) / "04-fluid-solid"
        os.chdir(here)

    config = here / "../dimensioning.yaml"
    mesh   = "mesh.msh"

    if not Path(mesh).exists():
        generate_charging_case(config, mesh)
    else:
        post_process(here)

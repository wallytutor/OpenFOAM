# -*- coding: utf-8 -*-
from ruamel.yaml import YAML
from heat_recovery.cowper_like import CowperLikeGeometry

yaml = YAML()
data = yaml.load(open("dimensioning.yaml"))

geom = CowperLikeGeometry(
    m_h = data["m_h"],
    D_h = data["D_h"],
    h_t = data["h_t"],

    # Use with care based on the above values:
    num_points_angular = 6,
    fluid_bl_tot = 0.003,
    fluid_bl_ext = 0.0001,
    fluid_bl_int = 0.0006,
    solid_bl_ext = 0.0050,
    solid_bl_int = 0.0001,
    rel_layer    = 0.15,
)

geom.create_model(saveas="mesh.msh", render=True)

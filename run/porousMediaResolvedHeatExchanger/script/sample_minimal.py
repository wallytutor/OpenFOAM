# -*- coding: utf-8 -*-
import os
import shutil
from pathlib import Path
import numpy as np
from numpy import sin, cos, pi

from porous_media.utils import get_blender
blender_exe = get_blender()

import trimesh
import pyvista as pv
import trimesh.interfaces.blender

import porous_media.geometry as pg

trimesh.interfaces.blender._blender_executable = blender_exe
trimesh.interfaces.blender.exists = True


def gyroid(X, Y, Z, **kwargs):
    """ Evaluate the gyroid implicit field on the input grid. """
    return sin(X) * cos(Y) + sin(Y) * cos(Z) + sin(Z) * cos(X)


def load_physical_mesh(filename, x_lims, y_lims, z_lims, nx, ny, nz):
    stl_trimesh = trimesh.load(filename)
    stl_physical = pg.map_voxel_to_physical(
        mesh   = stl_trimesh,
        x_lims = x_lims,
        y_lims = y_lims,
        z_lims = z_lims,
        nx     = n,
        ny     = n,
        nz     = n,
    )
    return stl_physical



if __name__ == "__main__":
    x_lims = (-pi, pi)
    y_lims = (-pi, pi)
    z_lims = (-pi, pi)
    n = 200
    box_frac = 0.8

    # Using 0.5001 instead of exactly 0.5 avoids saddle-point ambiguities in
    # marching cubes, ensuring a perfectly closed and watertight solid wall.
    surface = pg.FunctionalShapes(
        functional = gyroid,
        level      = 0.5001,
        thickness  = 1.0,
        x_lims     = x_lims,
        y_lims     = y_lims,
        z_lims     = z_lims,
        nx         = n,
        ny         = n,
        nz         = n,
    )

    filename = Path(__file__).parent / "surface.stl"
    surface.save_mesh(filename, show=False)

    stl_mesh = load_physical_mesh(
        filename = filename,
        x_lims   = x_lims,
        y_lims   = y_lims,
        z_lims   = z_lims,
        nx       = n,
        ny       = n,
        nz       = n,
    )

    subvolume = pg.get_subvolume(x_lims, y_lims, z_lims, box_frac)

    bodies = pg.get_porous_domain(stl_mesh, subvolume)

    pg.plot_domain(
        pore_bodies = bodies,
        # subvolume   = subvolume,
        # stl_mesh    = stl_mesh
    )

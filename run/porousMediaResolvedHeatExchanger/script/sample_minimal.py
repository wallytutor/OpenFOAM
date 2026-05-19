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


def map_voxel_to_physical(mesh, x_lims, y_lims, z_lims, nx, ny, nz):
    """ Scale mesh vertices from voxel index coordinates to physical coordinates.

    Marching cubes extracts vertices in index space [0, N-1], which needs
    to be mapped back to the physical x_lims, y_lims, z_lims range.
    """
    pts = mesh.vertices.copy()
    pts[:, 0] = x_lims[0] + pts[:, 0] * (x_lims[1] - x_lims[0]) / (nx - 1)
    pts[:, 1] = y_lims[0] + pts[:, 1] * (y_lims[1] - y_lims[0]) / (ny - 1)
    pts[:, 2] = z_lims[0] + pts[:, 2] * (z_lims[1] - z_lims[0]) / (nz - 1)
    mesh.vertices = pts
    return mesh


if __name__ == "__main__":
    x_lims = (-pi, pi)
    y_lims = (-pi, pi)
    z_lims = (-pi, pi)
    n = 200
    box_frac = 0.8

    surface = pg.FunctionalShapes(
        functional = gyroid,
        level      = 0.5,
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

    stl_trimesh = trimesh.load(filename)
    stl_physical = map_voxel_to_physical(
        mesh   = stl_trimesh,
        x_lims = x_lims,
        y_lims = y_lims,
        z_lims = z_lims,
        nx     = n,
        ny     = n,
        nz     = n,
    )


    box = pg.get_subvolume(x_lims, y_lims, z_lims, box_frac)
    pore_trimesh = box.difference(
        stl_physical, engine="blender", check_volume=False)

    print(f"Pore volume watertight? {pore_trimesh.is_watertight}")
    print(f"Pore volume exact value: {pore_trimesh.volume:.4f}")

    pore_pv = pv.wrap(pore_trimesh)
    bodies = pore_pv.split_bodies()

    print(f"Extracted pore space split into {len(bodies)} regions")
    for i, b in enumerate(bodies):
        print(f" - Region {i+1}: points={b.n_points}, cells={b.n_cells}, volume={b.volume:.4f}")

    pg.plot_domain(stl_physical, box, bodies)

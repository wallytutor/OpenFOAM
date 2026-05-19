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

    stl_physical = load_physical_mesh(
        filename = filename,
        x_lims   = x_lims,
        y_lims   = y_lims,
        z_lims   = z_lims,
        nx       = n,
        ny       = n,
        nz       = n,
    )

    box = pg.get_subvolume(x_lims, y_lims, z_lims, box_frac)

    pore_trimesh = box.difference(
        stl_physical, engine="blender", check_volume=False)

    # print(f"Solid wall mesh watertight? {stl_physical.is_watertight}")
    # print(f"Solid wall Euler number: {stl_physical.euler_number}")

    # print(f"Pore volume watertight? {pore_trimesh.is_watertight}")
    # print(f"Pore volume exact value: {pore_trimesh.volume:.4f}")

    # Use trimesh.split to split the pore mesh based on face adjacency,
    # which is topologically exact and unaffected by PyVista/VTK shared
    # vertex mergers.
    comps = pore_trimesh.split(only_watertight=False)

    min_vertices = 100
    bodies = []
    for c in comps:
        # Filter out tiny numerical boundary artifacts
        if len(c.vertices) < min_vertices:
            continue

        c.fix_normals()
        # Invert faces if volume is negative to ensure positive volume
        # in PyVista
        if c.volume < 0:
            c.invert()
        bodies.append(pv.wrap(c))

    # Sort largest components first
    bodies = sorted(bodies, key=lambda b: b.n_points, reverse=True)

    print(f"Extracted pore space split into {len(bodies)} regions")
    for i, b in enumerate(bodies):
        print(f" - Region {i+1}: points={b.n_points}, cells={b.n_cells}, volume={b.volume:.4f}")

    pg.plot_domain(stl_physical, box, bodies)

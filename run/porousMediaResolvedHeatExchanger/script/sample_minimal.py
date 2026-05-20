# -*- coding: utf-8 -*-
from pathlib import Path

import numpy as np
from numpy import sin, cos, pi

import porous_media.geometry as pg


def gyroid(X, Y, Z, **kwargs):
    """ Evaluate the gyroid implicit field on the input grid. """
    return sin(X) * cos(Y) + sin(Y) * cos(Z) + sin(Z) * cos(X)


if __name__ == "__main__":
    n        = 200
    x_lims   = (-pi, pi)
    y_lims   = (-pi, pi)
    z_lims   = (-pi, pi)
    box_frac = 0.9

    surface = pg.FunctionalShapes(
        functional = gyroid,
        level      = 0.500,
        thickness  = 5.0,
        x_lims     = x_lims,
        y_lims     = y_lims,
        z_lims     = z_lims,
        nx         = n,
        ny         = n,
        nz         = n,
    )

    filename = Path(__file__).parent / "surface.stl"
    surface.save_mesh(filename, show=False)

    # stl_mesh = pg.FunctionalShapes.from_stl(filename)
    # XXX the above only works if STL was created via FunctionalShapes!
    # Otherwise you may need to call `FunctionalShapes.map_voxel_to_physical`
    # on the mesh from STL to get the right voxel mapping!

    stl_mesh = surface.mesh
    x_lims, y_lims, z_lims = stl_mesh.bounding_box.bounds.T


    # subvolume = pg.get_subvolume(x_lims, y_lims, z_lims, box_frac)
    # bodies = pg.get_porous_domain(stl_mesh, subvolume)

    # pg.plot_domain(
    #     pore_bodies = bodies,
    #     # subvolume   = subvolume,
    #     # stl_mesh    = stl_mesh,
    # )

# -*- coding: utf-8 -*-

from pathlib import Path
from numpy import sin, cos, pi
from ruamel.yaml import YAML

import porous_media.geometry as pg

t = 0.002
""" Thickness of the solid walls [m]. """

L = 0.1
""" System overall dimensions [m]. """

k = 3
""" Number of wavelengths along each direction. """


def gyroid(X, Y, Z, **kwargs):
    """ Evaluate the gyroid implicit field on the input grid. """

    m = 2.0 * pi * k / L

    X_ = m * X
    Y_ = m * Y
    Z_ = m * Z

    return sin(X_) * cos(Y_) + sin(Y_) * cos(Z_) + sin(Z_) * cos(X_)


if __name__ == "__main__":
    box_frac = 0.9
    n        = 150

    x_lims   = (0.0, L)
    y_lims   = (0.0, L)
    z_lims   = (0.0, L)

    surface = pg.FunctionalShapes(
        functional = gyroid,
        level      = 0.500,
        thickness  = t,
        x_lims     = x_lims,
        y_lims     = y_lims,
        z_lims     = z_lims,
        nx         = n,
        ny         = n,
        nz         = n,
        centered   = False,
    )

    mesh = surface.mesh
    # volume = pg.FunctionalShapes.cube_subvolume(mesh, box_frac)
    # domain = pg.PorousDomainExtractor(mesh, volume)

    # mesh.metadata["type"] == "solid"

    # Save off-screen screenshot representation
    pg.plot_domain(
        bodies = [mesh],
        # bodies = domain.solid_bodies,
        # saveas = Path(__file__).parent / "geometry.png",
    )

    # Save extracted fluid and solid zones directly for snappyHexMesh
    # domain.save_project(
    #     project_dir  = Path(__file__).parent / "constant" / "triSurface",
    #     overwrite    = True,
    #     save_sources = False,
    # )

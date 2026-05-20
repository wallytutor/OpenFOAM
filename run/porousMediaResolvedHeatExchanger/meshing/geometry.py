# -*- coding: utf-8 -*-

from pathlib import Path
from numpy import sin, cos, pi
from ruamel.yaml import YAML

import porous_media.geometry as pg


def gyroid(X, Y, Z, **kwargs):
    """ Evaluate the gyroid implicit field on the input grid. """
    return sin(X) * cos(Y) + sin(Y) * cos(Z) + sin(Z) * cos(X)


if __name__ == "__main__":
    box_frac = 0.5
    n        = 100

    x_lims   = (-2 * pi, 2 * pi)
    y_lims   = (-2 * pi, 2 * pi)
    z_lims   = (-2 * pi, 2 * pi)

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

    mesh = surface.mesh
    volume = pg.FunctionalShapes.cube_subvolume(mesh, box_frac)
    domain = pg.PorousDomainExtractor(mesh, volume)

    # Save off-screen screenshot representation
    pg.plot_domain(
        bodies = domain.bodies,
        # saveas = Path(__file__).parent / "geometry.png",
    )

    # Save extracted fluid and solid zones directly for snappyHexMesh
    # domain.save_project(
    #     project_dir  = Path(__file__).parent / "constant" / "triSurface",
    #     overwrite    = True,
    #     save_sources = False,
    # )

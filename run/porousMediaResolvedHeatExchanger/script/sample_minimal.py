# -*- coding: utf-8 -*-
from pathlib import Path

import numpy as np
from numpy import sin, cos, pi

import porous_media.geometry as pg


def gyroid(X, Y, Z, **kwargs):
    """ Evaluate the gyroid implicit field on the input grid. """
    return sin(X) * cos(Y) + sin(Y) * cos(Z) + sin(Z) * cos(X)


if __name__ == "__main__":
    box_frac = 0.9
    inv_frac = 1.0 / box_frac
    n        = 200

    x_lims   = (-inv_frac * pi, inv_frac * pi)
    y_lims   = (-inv_frac * pi, inv_frac * pi)
    z_lims   = (-inv_frac * pi, inv_frac * pi)

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

    mesh   = pg.FunctionalShapes.from_stl(filename)
    volume = pg.FunctionalShapes.cube_subvolume(mesh, box_frac)
    domain = pg.PorousDomainExtractor(mesh, volume)

    pg.plot_domain(
        bodies = domain.fluid_bodies,
        volume = volume,
        parent = mesh,
        saveas = filename.with_suffix(".png"),
    )

    domain.save_project(
        project_dir  = filename.parent / "sample_minimal",
        overwrite    = True,
        save_sources = False,
    )

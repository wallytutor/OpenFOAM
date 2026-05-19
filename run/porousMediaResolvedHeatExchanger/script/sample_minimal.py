# -*- coding: utf-8 -*-
from numpy import sin, cos, pi
from pathlib import Path
from porous_media.geometry import FunctionalShapes


def gyroid(X, Y, Z, **kwargs):
    """ Evaluate the gyroid implicit field on the input grid. """
    return sin(X) * cos(Y) + sin(Y) * cos(Z) + sin(Z) * cos(X)


if __name__ == "__main__":
    surface = FunctionalShapes(
        functional = gyroid,
        level      = 0.5,
        thickness  = 1.0,
        x_lims     = (-pi, pi),
        y_lims     = (-pi, pi),
        z_lims     = (-pi, pi),
        nx         = 200,
        ny         = 200,
        nz         = 200,
    )

    filename = Path(__file__).parent / "surface.stl"
    surface.save_mesh(filename, show=True)

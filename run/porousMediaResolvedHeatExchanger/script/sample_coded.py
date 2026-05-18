# -*- coding: utf-8 -*-
import numpy as np
from pathlib import Path
from porous_media.geometry import FunctionalShapes

HERE = Path(__file__).parent
""" Path to the directory containing this script. """

n = 2
""" Number of porosity length scales in each dimension. """

l = 0.025
""" Porosity length scale (m). """

m = 2.0 * np.pi / l
""" Porosity frequency (rad/m). """


def implicit_surface(X, Y, Z, **kwargs):
    a = np.sin(m * X) * np.cos(m * Y)
    b = np.sin(m * Y) * np.cos(m * Z)
    c = np.sin(m * Z) * np.cos(m * X)
    return a + b + c


surface = FunctionalShapes(
    # Functional type (or callable):
    functional = implicit_surface,

    # Domain limits and resolution:
    x_lims        = (0.0, n * l),
    y_lims        = (0.0, n * l),
    z_lims        = (0.0, n * l),
    nx            = 50,
    ny            = 50,
    nz            = 50,

    # Functional parameters:
    # rugosity_ampl = 0.005,
    # stripes_ampl  = 0.001,
    # stripes_freq  = 25,
    # level         = 15,
)

filename = HERE / "geometry.stl"
surface.save_mesh(filename, show=True)

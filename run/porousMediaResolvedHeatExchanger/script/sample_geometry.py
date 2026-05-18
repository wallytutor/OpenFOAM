# -*- coding: utf-8 -*-
import numpy as np
from pathlib import Path
from porous_media.geometry import FunctionalShapes

HERE = Path(__file__).parent
""" Path to the directory containing this script. """

surface = FunctionalShapes(
    # Functional type (or callable):
    functional = "schwarz",

    # Domain limits and resolution:
    x_lims        = (-3*np.pi, 3*np.pi),
    y_lims        = (-3*np.pi, 3*np.pi),
    z_lims        = (       0, 3*np.pi),
    nx            = 50,
    ny            = 50,
    nz            = 50,

    # Functional parameters:
    rugosity_ampl = 0.0,
    stripes_ampl  = 0.0,
    stripes_freq  = 0.0,
    level         = 0.0,

    # Functional-specific parameters:
    isocontour    = 0.0,
)

filename = HERE / "surface.stl"
surface.save_mesh(filename, show=True)
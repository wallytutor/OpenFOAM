# -*- coding: utf-8 -*-
import numpy as np
from pathlib import Path
from typing import Any
from numpy.typing import NDArray
from porous_media.geometry import FunctionalShapes
from porous_media.geometry import decorate_surface

# ---------------------------------------------------------------------
# User-defined implicit surfaces
# ---------------------------------------------------------------------

@decorate_surface
def gyroid(
        X: NDArray[Any],
        Y: NDArray[Any],
        Z: NDArray[Any],
        **kwargs: Any,
    ) -> NDArray[Any]:
    """ Evaluate the gyroid implicit field on the input grid.

    https://en.wikipedia.org/wiki/Gyroid
    """
    mx = kwargs.get("mx", 1)
    my = kwargs.get("my", 1)
    mz = kwargs.get("mz", 1)

    a = np.sin(mx * X) * np.cos(my * Y)
    b = np.sin(my * Y) * np.cos(mz * Z)
    c = np.sin(mz * Z) * np.cos(mx * X)

    return a + b + c


@decorate_surface
def schwarz_p(
        X: NDArray[Any],
        Y: NDArray[Any],
        Z: NDArray[Any],
        **kwargs: Any,
    ) -> NDArray[Any]:
    """ Evaluate the Schwarz-P implicit field on the input grid.

    https://en.wikipedia.org/wiki/Schwarz_minimal_surface
    """
    mx = kwargs.get("mx", 1)
    my = kwargs.get("my", 1)
    mz = kwargs.get("mz", 1)
    contour = kwargs.get("contour", 0.0)

    return np.cos(mx * X) + np.cos(my * Y) + np.cos(mz * Z) - contour

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

HERE = Path(__file__).parent

SMOOTH = {}

ROUGH = {
    "relative_roughness": 0.005,
    "wave_amplitude": 0.005,
    "wave_number": 0.1
}

WAVENUMBERS = {
    "mx": 3,
    "my": 3,
    "mz": 3,
}

# ---------------------------------------------------------------------
# Workflow function
# ---------------------------------------------------------------------

def generate_surface(func, config, saveas, **kwargs):
    """ Generate and save surface mesh for the given functional. """
    n = kwargs.get("n", 200)

    surface = FunctionalShapes(
        functional = func,
        level      = kwargs.get("level", 0.5),
        thickness  = kwargs.get("thickness", 1.0),
        x_lims     = (-np.pi, np.pi),
        y_lims     = (-np.pi, np.pi),
        z_lims     = (-np.pi, np.pi),
        nx         = n,
        ny         = n,
        nz         = n,
        **config,
        **WAVENUMBERS
    )

    filename = HERE / saveas
    surface.save_mesh(filename, show=True)


if __name__ == "__main__":
    generate_surface(gyroid,    SMOOTH, "gyroid_smooth.stl")
    generate_surface(gyroid,    ROUGH,  "gyroid_rough.stl")
    generate_surface(schwarz_p, SMOOTH, "schwarz_p_smooth.stl")
    generate_surface(schwarz_p, ROUGH,  "schwarz_p_rough.stl")

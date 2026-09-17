# -*- coding: utf-8 -*-
import numpy as np
from typing import Any
from numpy.typing import NDArray
from porous_media.geometry import decorate_surface


@decorate_surface
def implicit_surface(
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

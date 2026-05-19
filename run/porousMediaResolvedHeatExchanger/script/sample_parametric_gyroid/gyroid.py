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
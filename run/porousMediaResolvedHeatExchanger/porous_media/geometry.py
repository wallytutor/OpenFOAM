# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from skimage.measure import marching_cubes
from typing import Callable
from trimesh import Trimesh


def extrema(a: NDArray) -> tuple[np.number, np.number]:
    """ Calculate the minimum and maximum values of an array. """
    return a.min(), a.max()


def meshgrid_3d(
        *,
        x_lims: tuple[float, float],
        y_lims: tuple[float, float],
        z_lims: tuple[float, float],
        nx: int,
        ny: int,
        nz: int,
        sx: float = 1.0,
        sy: float = 1.0,
        sz: float = 1.0
    ) -> tuple[tuple[NDArray, NDArray, NDArray], NDArray]:
    """ Generate a 3D mesh grid within specified limits.

    Parameters
    ----------
    x_lims : tuple[float, float]
        The limits for the x-axis (min, max).
    y_lims : tuple[float, float]
        The limits for the y-axis (min, max).
    z_lims : tuple[float, float]
        The limits for the z-axis (min, max).
    nx : int
        The number of points along the x-axis.
    ny : int
        The number of points along the y-axis.
    nz : int
        The number of points along the z-axis.
    sx : float, optional
        Scaling factor for the x-axis (default is 1.0).
    sy : float, optional
        Scaling factor for the y-axis (default is 1.0).
    sz : float, optional
        Scaling factor for the z-axis (default is 1.0).
    """
    x = sx * np.linspace(*x_lims, nx)
    y = sy * np.linspace(*y_lims, ny)
    z = sz * np.linspace(*z_lims, nz)

    grid = np.meshgrid(x, y, z)
    scale = np.array([sx, sy, sz])

    return grid, scale


class FunctionalShapes:
    def __init__(self,
            X: NDArray[np.number],
            Y: NDArray[np.number],
            Z: NDArray[np.number],
            *,
            functional: str | Callable,
            rugosity_ampl: float =  0.0,
            stripes_ampl: float = 0.0,
            stripes_freq: float = 1.0,
            **kwargs
        ):
        if isinstance(functional, str):
            f = getattr(self, functional)
        else:
            f = functional

        surface = f(X, Y, Z, **kwargs)

        if rugosity_ampl > 0.0:
            surface += rugosity_ampl * np.random.normal(size=surface.shape)

        if stripes_ampl > 0.0 and stripes_freq > 0.0:
            surface += stripes_ampl * np.sin(stripes_freq * Z)

        vertices, faces = self._extract_features(surface, **kwargs)

        self._vertices = vertices
        self._faces    = faces
        self._mesh     = Trimesh(vertices=vertices, faces=faces)

    @staticmethod
    def _extract_features(surface, level=0, **kwargs):
        vmin, vmax = extrema(surface)

        if not (vmin <= level <= vmax):
            raise ValueError(f"Level {level} outside [{vmin}; {vmax}]")

        vertices, faces, _, _ = marching_cubes(surface, level=level)
        faces = faces.astype(np.int32)

        if len(vertices) == 0 or len(faces) == 0:
            raise RuntimeError("No surface generated!")

        return vertices, faces

    @classmethod
    def gyroid(cls, X, Y, Z, **kwargs):
        a = np.sin(X) * np.cos(Y)
        b = np.sin(Y) * np.cos(Z)
        c = np.sin(Z) * np.cos(X)
        return a + b + c

    @classmethod
    def schwarz(cls, X, Y, Z, isocontour=0.0, **kwargs):
        return np.cos(X) + np.cos(Y) + np.cos(Z) - isocontour

    @property
    def vertices(self):
        return self._vertices

    @property
    def faces(self):
        return self._faces

    @property
    def mesh(self):
        return self._mesh

    def save_mesh(self, filename: str):
        self._mesh.export(file_obj=filename)

    def plot(self, scale, **kwargs):
        vx = self._vertices[:, 0] * scale[0]
        vy = self._vertices[:, 1] * scale[1]
        vz = self._vertices[:, 2] * scale[2]

        opts = dict(
            cmap  = kwargs.get("cmap", "gnuplot"),
            lw    = kwargs.get("lw", 0.1),
            alpha = kwargs.get("alpha", 1.0)
        )

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(projection="3d")
        ax.plot_trisurf(vx, vy, vz, triangles=self._faces, **opts)
        ax.auto_scale_xyz(extrema(vx), extrema(vy), extrema(vz))
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        fig.tight_layout()

        return fig, ax

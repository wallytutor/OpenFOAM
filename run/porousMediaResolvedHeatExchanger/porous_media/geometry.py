# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from numpy.typing import NDArray
from skimage.measure import marching_cubes
from typing import Any, Protocol
from matplotlib.figure import Figure
from trimesh import Trimesh


class FunctionalType(Protocol):
    """ Protocol for a functional that generates a 3D surface. """
    def __call__(
        self,
        X: NDArray[Any],
        Y: NDArray[Any],
        Z: NDArray[Any],
        **kwargs: Any
    ) -> NDArray[Any]:
        ...


class FunctionalShapes:
    """ Generate triangulated 3D isosurfaces from analytic scalar fields.

    The class evaluates a scalar field on a regular grid and extracts an
    isosurface with marching cubes. The resulting geometry is stored as raw
    arrays (`vertices`, `faces`) and as a `trimesh.Trimesh` object (`mesh`).

    Parameters
    ----------
    functional : str | FunctionalType
        Name of a built-in functional (for example ``"gyroid"`` or
        ``"schwarz"``) or a callable with signature
        ``f(X, Y, Z, **kwargs) -> NDArray``.
    x_lims : tuple[float, float]
        Domain limits for x as ``(xmin, xmax)``.
    y_lims : tuple[float, float]
        Domain limits for y as ``(ymin, ymax)``.
    z_lims : tuple[float, float]
        Domain limits for z as ``(zmin, zmax)``.
    nx : int, default=100
        Number of grid points along x.
    ny : int, default=100
        Number of grid points along y.
    nz : int, default=100
        Number of grid points along z.
    rugosity_ampl : float, default=0.0
        Amplitude of additive Gaussian noise applied to the sampled field.
    stripes_ampl : float, default=0.0
        Amplitude of sinusoidal modulation added along z.
    stripes_freq : float, default=1.0
        Frequency used for sinusoidal stripe modulation.
    level : float, default=0.0
        Isovalue extracted from the sampled scalar field.
    **kwargs : Any
        Extra keyword arguments forwarded to the selected functional.

    Raises
    ------
    ValueError
        If the functional is invalid or if ``level`` lies outside the
        sampled field range.
    RuntimeError
        If marching cubes cannot extract any surface.
    """
    def __init__(self, *,
            functional: str | FunctionalType,
            x_lims: tuple[float, float],
            y_lims: tuple[float, float],
            z_lims: tuple[float, float],
            nx: int = 100,
            ny: int = 100,
            nz: int = 100,
            rugosity_ampl: float =  0.0,
            stripes_ampl: float = 0.0,
            stripes_freq: float = 1.0,
            level: float = 0.0,
            **kwargs: Any,
        ) -> None:
        """Initialize the surface by sampling a field and extracting an isosurface."""
        #region: generate grid
        (X, Y, Z), _ = self.meshgrid_3d(x_lims=x_lims, y_lims=y_lims,
                                        z_lims=z_lims, nx=nx, ny=ny, nz=nz)
        #endregion: generate grid

        #region: resolve functional
        if callable(functional):
            f = functional
        elif isinstance(functional, str):
            try:
                f = getattr(self, functional)
            except AttributeError:
                raise ValueError(f"Functional '{functional}' not found")
        else:
            raise ValueError("Functional must be a string or a callable")
        #endregion: resolve functional

        #region: generate surface
        surface = f(X, Y, Z, **kwargs)

        if rugosity_ampl > 0.0:
            surface += rugosity_ampl * np.random.normal(size=surface.shape)

        if stripes_ampl > 0.0 and stripes_freq > 0.0:
            surface += stripes_ampl * np.sin(stripes_freq * Z)
        #endregion: generate surface

        #region: extract features
        vmin, vmax = surface.min(), surface.max()

        if not (vmin <= level <= vmax):
            raise ValueError(f"Level {level} outside [{vmin}; {vmax}]")

        vertices, faces, _, _ = marching_cubes(surface, level=level)
        faces = faces.astype(np.int32)

        if len(vertices) == 0 or len(faces) == 0:
            raise RuntimeError("No surface generated!")
        #endregion: extract features

        self._vertices = vertices
        self._faces    = faces
        self._mesh     = Trimesh(vertices=vertices, faces=faces)

    @classmethod
    def gyroid(
            cls,
            X: NDArray[Any],
            Y: NDArray[Any],
            Z: NDArray[Any],
            **kwargs: Any,
        ) -> NDArray[Any]:
        """ Evaluate the gyroid implicit field on the input grid. """
        a = np.sin(X) * np.cos(Y)
        b = np.sin(Y) * np.cos(Z)
        c = np.sin(Z) * np.cos(X)
        return a + b + c

    @classmethod
    def schwarz(
            cls,
            X: NDArray[Any],
            Y: NDArray[Any],
            Z: NDArray[Any],
            isocontour: float = 0.0,
            **kwargs: Any,
        ) -> NDArray[Any]:
        """ Evaluate the Schwarz-P implicit field on the input grid. """
        return np.cos(X) + np.cos(Y) + np.cos(Z) - isocontour

    @staticmethod
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

    @property
    def vertices(self) -> NDArray[Any]:
        """ Return mesh vertices with shape ``(n_vertices, 3)``. """
        return self._vertices

    @property
    def faces(self) -> NDArray[Any]:
        """ Return triangular face indices with shape ``(n_faces, 3)``. """
        return self._faces

    @property
    def mesh(self) -> Trimesh:
        """ Return the `trimesh.Trimesh` view of the generated surface. """
        return self._mesh

    def save_mesh(self, filename: str) -> None:
        """ Export the current mesh to a file path supported by trimesh. """
        self._mesh.export(file_obj=filename)

    def plot_mesh(self, filename: str) -> None:
        """ Load a mesh file with PyVista display STL file. """
        mesh = pv.read(filename)
        mesh.plot()

    def render_pyplot(self, **kwargs: Any) -> tuple[Figure, Any]:
        """ Render the generated surface with Matplotlib. """
        vx = self._vertices[:, 0]
        vy = self._vertices[:, 1]
        vz = self._vertices[:, 2]

        dx = vx.min(), vx.max()
        dy = vy.min(), vy.max()
        dz = vz.min(), vz.max()

        opts = dict(
            cmap  = kwargs.get("cmap", "gnuplot"),
            lw    = kwargs.get("lw", 0.1),
            alpha = kwargs.get("alpha", 1.0)
        )

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(projection="3d")
        ax.plot_trisurf(vx, vy, vz, triangles=self._faces, **opts)
        ax.auto_scale_xyz(dx, dy, dz)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        fig.tight_layout()

        return fig, ax

# -*- coding: utf-8 -*-

import argparse
import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from numpy.typing import NDArray
from ruamel.yaml import YAML
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
        Name of a built-in functional (``"gyroid"`` or ``"schwarz"``),
        a path to a Python module file (any string containing a ``.``
        before a recognised extension such as ``.py``) that exposes
        ``implicit_surface(X, Y, Z, **kwargs) -> NDArray``, or a
        plain callable with the same signature.
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
        #region: generate grid
        (X, Y, Z), _ = self.meshgrid_3d(x_lims=x_lims, y_lims=y_lims,
                                        z_lims=z_lims, nx=nx, ny=ny, nz=nz)
        #endregion: generate grid

        #region: generate surface
        f = self._resolve_functional(functional)
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

    def _resolve_functional(self,
            functional: str | FunctionalType
        ) -> FunctionalType:
        """ Resolve a functional to be used for mesh generation. """
        if callable(functional):
            return functional

        if isinstance(functional, str):
            # treat as a module if any suffix is present
            if (path := Path(functional)).suffix:
                return self._load_module_functional(path)
            else:
                try:
                    return getattr(self, functional)
                except AttributeError:
                    raise ValueError(f"Functional '{functional}' not found")

        raise ValueError("Functional must be a string or a callable")

    @staticmethod
    def _load_module_functional(path: Path) -> FunctionalType:
        """ Import ``implicit_surface``  from path to module.

        Parameters
        ----------
        path : Path
            Absolute or relative path to a Python source file that defines
            ``implicit_surface(X, Y, Z, **kwargs) -> NDArray``.

        Raises
        ------
        FileNotFoundError
            If *path* does not point to an existing file.
        AttributeError
            If the loaded module does not expose ``implicit_surface``.
        """
        resolved = path.resolve()

        if not resolved.is_file():
            raise FileNotFoundError(f"Functional module not found: {resolved}")

        spec = importlib.util.spec_from_file_location(resolved.stem, resolved)

        if spec is None or spec.loader is None:
            raise ImportError(f"Cannot load module from: {resolved}")

        module: ModuleType = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)  # type: ignore[union-attr]

        if not hasattr(module, "implicit_surface"):
            raise AttributeError(
                f"Module '{resolved}' does not define 'implicit_surface'"
            )

        return module.implicit_surface  # type: ignore[return-value]

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


def main() -> int:
    #region: build parser
    parser = argparse.ArgumentParser(
        prog="porous-media-geometry",
        description="Generate porous media geometry from YAML configuration.",
    )

    parser.add_argument(
        "config",
        type=Path,
        help="Path to YAML file with FunctionalShapes parameters.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting existing output files (default: false).",
    )

    args = parser.parse_args()
    #endregion: build parser

    #region: load config
    yaml = YAML(typ="safe")

    with args.config.open("r", encoding="utf-8") as fobj:
        config = yaml.load(fobj)

    if config is None:
        raise ValueError(f"Empty YAML file: {args.config}")

    if not isinstance(config, dict):
        raise ValueError("Top-level YAML content must be a mapping")

    if not isinstance(params := config["functional_shapes"], dict):
        raise ValueError("'functional_shapes' must be a mapping")
    #endregion: load config

    #region: validate output
    plot_mesh = bool(config.get("plot_mesh", True))
    output_path = Path(config.get("output", "surface.stl"))

    if output_path.exists() and not args.overwrite:
        raise FileExistsError(f"Output file already exists: {output_path}")

    if output_path.suffix.lower() != ".stl":
        raise ValueError(f"Unsupported output format: {output_path.suffix}")
    #endregion: validate output

    surface = FunctionalShapes(**params)
    surface.save_mesh(str(output_path))

    if plot_mesh:
        surface.plot_mesh(str(output_path))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

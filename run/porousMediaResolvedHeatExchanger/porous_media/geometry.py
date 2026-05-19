# -*- coding: utf-8 -*-

import argparse
import functools
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
    level : float, default=0.0
        Fractional level for isosurface extraction, where 0.0 corresponds
        to the minimum value of the sampled field and 1.0 to the maximum.
        For example, a value of 0.5 will extract the isosurface at the
        midpoint between the minimum and maximum field values.
    **kwargs : Any
        Extra keyword arguments forwarded to the selected functional.

    Raises
    ------
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
            level: float = 0.0,
            **kwargs: Any,
        ) -> None:
        f = self._resolve_functional(functional)

        # XXX keep as a function (we may need hexagonal grids later!)
        X, Y, Z = self.cube_meshgrid(x_lims, y_lims, z_lims, nx, ny, nz)

        surface = f(X, Y, Z, **kwargs)
        vertices, faces = self._extract_features(surface, level)

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
            if not (module := Path(functional)).exists():
                raise FileNotFoundError(f"No such module {functional}")

            return load_implicit_surface(module)

        raise ValueError("Functional must be a string or a callable")

    def _extract_features(self,
            surface: NDArray[Any],
            fractional_level: float
        ) -> tuple[NDArray, NDArray]:
        """ Extract features from a scalar field using marching cubes. """
        vmin, vmax = surface.min(), surface.max()

        level = vmin + fractional_level * (vmax - vmin)

        vertices, faces, _, _ = marching_cubes(surface, level=level)

        faces = faces.astype(np.int32)

        if len(vertices) == 0 or len(faces) == 0:
            raise RuntimeError("No surface generated!")

        return vertices, faces

    @staticmethod
    def cube_meshgrid(
            x_lims: tuple[float, float],
            y_lims: tuple[float, float],
            z_lims: tuple[float, float],
            nx: int,
            ny: int,
            nz: int
        ) -> tuple[NDArray, NDArray, NDArray]:
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
        """
        x = np.linspace(*x_lims, nx)
        y = np.linspace(*y_lims, ny)
        z = np.linspace(*z_lims, nz)
        return np.meshgrid(x, y, z)

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

    def save_mesh(self, filename: str | Path, show: bool = False) -> None:
        """ Export the current mesh to a file path supported by trimesh. """
        self._mesh.export(file_obj=filename)

        if show:
            self.plot_mesh(filename)

    def plot_mesh(self, filename: str | Path) -> None:
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


def noisy_surface(
        surface: NDArray[np.number],
        rugosity_ampl: float
    ) -> NDArray[np.number]:
    """ Apply additive Gaussian noise to a surface.

    Parameters
    ----------
    surface : NDArray[np.number]
        The input surface array.
    rugosity_ampl : float
        Amplitude of the Gaussian noise.

    Returns
    -------
    NDArray[np.number]
        The surface with added noise.
    """
    if rugosity_ampl <= 0.0:
        return surface

    noise = np.random.normal(size=surface.shape)
    return surface + rugosity_ampl * noise


def wavy_surface(
        surface: NDArray[np.number],
        Z: NDArray[np.number],
        amplitude: float,
        wave_number: float
    ) -> NDArray[np.number]:
    """ Apply sinusoidal modulation along the z-axis to a surface.

    Parameters
    ----------
    surface : NDArray[np.number]
        The input surface array.
    Z : NDArray[np.number]
        The z-coordinates corresponding to the surface points.
    amplitude : float
        Amplitude of the sinusoidal modulation.
    wave_number : float
        Wave number of the sinusoidal modulation.

    Returns
    -------
    NDArray[np.number]
        The surface with added sinusoidal modulation.
    """
    if amplitude <= 0.0 or wave_number <= 0.0:
        return surface

    return surface + amplitude * np.sin(2.0 * np.pi / wave_number * Z)


def decorate_surface(func: FunctionalType) -> FunctionalType:
    """ Decorator to apply noise and/or sinusoidal modulation to a surface. """
    @functools.wraps(func)
    def decorated(
            X: NDArray[Any],
            Y: NDArray[Any],
            Z: NDArray[Any],
            **kwargs: Any
        ) -> NDArray[Any]:
        A = kwargs.get("relative_roughness", 0.0)
        B = kwargs.get("wave_amplitude", 0.0)
        f = kwargs.get("wave_number", 0.0)

        surface = func(X, Y, Z, **kwargs)
        vmin, vmax = surface.min(), surface.max()
        # print(f"Functional range: [{vmin:.3f}, {vmax:.3f}]")

        surface = noisy_surface(surface, A * abs(vmax - vmin))
        vmin, vmax = surface.min(), surface.max()
        # print(f"Functional range: [{vmin:.3f}, {vmax:.3f}]")

        surface = wavy_surface(surface, Z, B * abs(vmax - vmin), f)
        # vmin, vmax = surface.min(), surface.max()
        # print(f"Functional range: [{vmin:.3f}, {vmax:.3f}]")

        return surface

    return decorated


def load_implicit_surface(path: Path) -> FunctionalType:
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
    surface.save_mesh(str(output_path), show=plot_mesh)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

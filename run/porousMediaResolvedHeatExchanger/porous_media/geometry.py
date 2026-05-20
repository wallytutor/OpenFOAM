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
import trimesh

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
    thickness : float, default=0.0
        Optional shell thickness. If positive, the extracted surface is
        transformed into a thickened shell using vertex-normal offsets.
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
            thickness: float = 0.0,
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

        if thickness > 0.0:
            self.thicken(thickness=thickness, inplace=True)

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

    def thicken(self,
            thickness: float,
            *,
            centered: bool = True,
            inplace: bool = False,
        ) -> Trimesh:
        """Create a thickened shell mesh from the current surface.

        Parameters
        ----------
        thickness : float
            Total shell thickness.
        centered : bool, default=True
            If true, offset half the thickness on each side of the current
            surface. If false, keep the current surface as the outer side and
            offset only inward.
        inplace : bool, default=False
            If true, replace the current mesh/vertices/faces with the
            thickened shell.

        Returns
        -------
        Trimesh
            Thickened shell mesh.
        """
        thickened = thicken_surface_mesh(
            self._mesh,
            thickness=thickness,
            centered=centered,
        )

        if inplace:
            self._mesh = thickened
            self._vertices = thickened.vertices
            self._faces = thickened.faces

        return thickened

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


def _boundary_edges(faces: NDArray[Any]) -> NDArray[Any]:
    """ Return edges that belong to exactly one triangle. """
    edges = np.vstack((
        faces[:, [0, 1]],
        faces[:, [1, 2]],
        faces[:, [2, 0]],
    ))

    sorted_edges = np.sort(edges, axis=1)
    unique_edges, counts = np.unique(sorted_edges, axis=0, return_counts=True)
    return unique_edges[counts == 1]


def thicken_surface_mesh(
        mesh: Trimesh,
        *,
        thickness: float,
        centered: bool = True,
    ) -> Trimesh:
    """ Create a shell with finite thickness from a triangulated surface.

    The algorithm offsets vertices along vertex normals to form two surfaces
    and connects boundary edges (if present) with side triangles.

    Parameters
    ----------
    mesh : Trimesh
        Input surface mesh.
    thickness : float
        Total shell thickness (must be positive).
    centered : bool, default=True
        If true, offsets equally to both sides of the original surface.
        If false, keeps original surface as the outer side and offsets inward.

    Returns
    -------
    Trimesh
        Thickened shell mesh.
    """
    if not np.isfinite(thickness):
        raise ValueError(f"thickness must be finite, got {thickness}")

    if thickness <= 0.0:
        raise ValueError(f"thickness must be positive, got {thickness}")

    vertices = np.asarray(mesh.vertices)
    faces = np.asarray(mesh.faces, dtype=np.int32)
    normals = np.asarray(mesh.vertex_normals)

    if len(vertices) == 0 or len(faces) == 0:
        raise ValueError("Cannot thicken an empty mesh")

    if centered:
        outer_shift = 0.5 * thickness
        inner_shift = -0.5 * thickness
    else:
        outer_shift = 0.0
        inner_shift = -thickness

    outer_vertices = vertices + outer_shift * normals
    inner_vertices = vertices + inner_shift * normals

    n_vertices = len(vertices)
    outer_faces = faces.copy()
    inner_faces = faces[:, ::-1] + n_vertices

    boundary = _boundary_edges(faces)
    if len(boundary) > 0:
        a = boundary[:, 0]
        b = boundary[:, 1]

        side_1 = np.column_stack((a, b, b + n_vertices))
        side_2 = np.column_stack((a, b + n_vertices, a + n_vertices))
        side_faces = np.vstack((side_1, side_2)).astype(np.int32)
        shell_faces = np.vstack((outer_faces, inner_faces, side_faces))
    else:
        shell_faces = np.vstack((outer_faces, inner_faces))

    shell_vertices = np.vstack((outer_vertices, inner_vertices))
    return Trimesh(vertices=shell_vertices, faces=shell_faces, process=True)


def map_voxel_to_physical(mesh, x_lims, y_lims, z_lims, nx, ny, nz):
    """ Scale mesh vertices from voxel space to physical coordinates.

    Marching cubes extracts vertices in index space [0, N-1], which needs
    to be mapped back to the physical x_lims, y_lims, z_lims range.
    Note: We must account for the 'xy' indexing of np.meshgrid in
    FunctionalShapes.cube_meshgrid, which swaps the X and Y axes.
    """
    delta_x = (x_lims[1] - x_lims[0]) / (nx - 1)
    delta_y = (y_lims[1] - y_lims[0]) / (ny - 1)
    delta_z = (z_lims[1] - z_lims[0]) / (nz - 1)

    pts = mesh.vertices.copy()
    pts_physical = np.zeros_like(pts)
    pts_physical[:, 0] = x_lims[0] + pts[:, 1] * delta_x
    pts_physical[:, 1] = y_lims[0] + pts[:, 0] * delta_y
    pts_physical[:, 2] = z_lims[0] + pts[:, 2] * delta_z
    mesh.vertices = pts_physical

    return mesh


def get_subvolume(x_lims, y_lims, z_lims, box_frac=0.8):
    """ Create a hexahedral/cubic region (trimesh box) that is a fraction of the bounding box.

    Parameters
    ----------
    x_lims : tuple[float, float]
        Original domain limits for x.
    y_lims : tuple[float, float]
        Original domain limits for y.
    z_lims : tuple[float, float]
        Original domain limits for z.
    box_frac : float, default=0.8
        Fraction of the original bounding box dimensions to scale by.

    Returns
    -------
    trimesh.Trimesh
        Sub-volume box centered on the original bounding box.
    """
    cx = (x_lims[0] + x_lims[1]) / 2.0
    cy = (y_lims[0] + y_lims[1]) / 2.0
    cz = (z_lims[0] + z_lims[1]) / 2.0

    dx = (x_lims[1] - x_lims[0]) * box_frac
    dy = (y_lims[1] - y_lims[0]) * box_frac
    dz = (z_lims[1] - z_lims[0]) * box_frac

    box = trimesh.creation.box(extents=(dx, dy, dz))
    box.apply_translation([cx, cy, cz])

    return box


def _filter_bodies(
        bodies: list[trimesh.Trimesh],
        min_vertices: int = 100,
    ) -> list[trimesh.Trimesh]:

    result = []

    for c in bodies:
        if len(c.vertices) < min_vertices:
            continue

        c.fix_normals()
        # Invert faces if volume is negative to ensure positive
        # volume in PyVista
        if c.volume < 0:
            c.invert()

        result.append(c)

    # Sort largest components first
    return result


def get_porous_domain(
        stl_physical: trimesh.Trimesh,
        stl_domain: trimesh.Trimesh,
        min_vertices: int = 100,
    ) -> list[trimesh.Trimesh]:

    # 1. Fluid domain: subtract the solid wall from the sub-volume
    pore_trimesh = stl_domain.difference(
        stl_physical, engine="blender", check_volume=False)
    pore_bodies = pore_trimesh.split(only_watertight=False)
    pore_bodies = _filter_bodies(pore_bodies, min_vertices=min_vertices)
    for c in pore_bodies:
        c.metadata["type"] = "fluid"

    # Sort largest fluid components first
    pore_bodies = sorted(pore_bodies, key=lambda b: len(b.vertices), reverse=True)

    # 2. Solid domain (the gap): intersect the solid wall with the sub-volume
    solid_trimesh = stl_domain.intersection(
        stl_physical, engine="blender", check_volume=False)
    solid_bodies = solid_trimesh.split(only_watertight=False)
    solid_bodies = _filter_bodies(solid_bodies, min_vertices=min_vertices)
    for c in solid_bodies:
        c.metadata["type"] = "solid"

    # Sort largest solid components first
    solid_bodies = sorted(solid_bodies, key=lambda b: len(b.vertices), reverse=True)

    return pore_bodies + solid_bodies


def plot_domain(
        pore_bodies,
        subvolume: trimesh.Trimesh | None = None,
        stl_mesh: trimesh.Trimesh | None = None,
    ) -> None:
    """ Display the extracted domain components.

    Parameters
    ----------
    pore_bodies : list[trimesh.Trimesh]
        List of extracted domain components (fluid and solid).
    subvolume : trimesh.Trimesh
        Sub-volume box mesh.
    stl_mesh : trimesh.Trimesh
        Original solid domain mesh.
    """
    print("Rendering final zones in PyVista...")
    plotter = pv.Plotter(title="Watertight Pore Volume Zones")

    # Add thin reference model of the solid wall shell
    if stl_mesh is not None:
        solid_pv = pv.wrap(stl_mesh)
        plotter.add_mesh(
            mesh        = solid_pv,
            color       = "#888888",
            opacity     = 0.12,
            show_edges  = False,
            label       = "Solid Wall (Reference)"
        )

    # Add wireframe bounding box
    if subvolume is not None:
        box_pv = pv.wrap(subvolume)
        plotter.add_mesh(
            mesh        = box_pv,
            color       = "#ffffff",
            style       = "wireframe",
            line_width  = 2,
            label       = "Sub-Volume Bounds"
        )

    # Vibrant color palettes
    fluid_colors = ["#3A86FF", "#FF006E", "#8338EC", "#06D6A0", "#FB5607"]
    solid_colors = ["#FFBE0B", "#E76F51", "#2A9D8F", "#F4A261"]

    fluid_idx = 0
    solid_idx = 0

    for body in pore_bodies:
        body_type = body.metadata.get("type", "fluid")
        if body_type == "fluid":
            color = fluid_colors[fluid_idx % len(fluid_colors)]
            label = f"Pore Zone {fluid_idx + 1} (Vol: {body.volume:.4f})"
            opacity = 0.85
            fluid_idx += 1
        else:
            color = solid_colors[solid_idx % len(solid_colors)]
            label = f"Solid Zone {solid_idx + 1} (Vol: {body.volume:.4f})"
            opacity = 0.50  # Semi-translucent to see interpenetrating networks
            solid_idx += 1

        plotter.add_mesh(
            mesh        = body,
            color       = color,
            opacity     = opacity,
            show_edges  = False,
            label       = label
        )

    plotter.add_legend(face="circle", bcolor=None)
    plotter.show_axes()
    plotter.show()



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

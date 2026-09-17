# -*- coding: utf-8 -*-

from IPython import display
import argparse
import os
import functools
import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Any, Self

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

# This must come before trimesh!
from .utils import get_blender
blender_exe = get_blender()

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
    """ Generate 3D isosurfaces from scalar fields.

    The class evaluates a scalar field on a regular grid and
    extracts an isosurface with marching cubes.

    Parameters
    ----------
    functional : str | FunctionalType
        Name of a built-in functional ("gyroid" or "schwarz"),
        a path to a Python module file, or a callable.
    x_lims : tuple[float, float]
        Domain limits for x as (xmin, xmax).
    y_lims : tuple[float, float]
        Domain limits for y as (ymin, ymax).
    z_lims : tuple[float, float]
        Domain limits for z as (zmin, zmax).
    nx : int, default=100
        Number of grid points along x.
    ny : int, default=100
        Number of grid points along y.
    nz : int, default=100
        Number of grid points along z.
    level : float, default=0.0
        Fractional level for isosurface extraction.
    thickness : float, default=0.0
        Optional shell thickness.
    centered : bool, default=False
        Whether to center the sub-volume on the original bounding box.
    **kwargs : Any
        Extra keyword arguments.
    """

    __slots__ = (
        "_vertices",
        "_faces",
        "_mesh",
    )

    def __new__(cls, *args: Any, **kwargs: Any) -> Self:
        obj = super().__new__(cls)
        obj._vertices = None
        obj._faces    = None
        obj._mesh     = None
        return obj

    def __init__(
            self,
            *,
            functional: str | FunctionalType,
            x_lims: tuple[float, float],
            y_lims: tuple[float, float],
            z_lims: tuple[float, float],
            nx: int = 100,
            ny: int = 100,
            nz: int = 100,
            level: float = 0.0,
            thickness: float = 0.0,
            centered: bool = False,
            **kwargs: Any,
        ) -> None:
        f = self._resolve_functional(functional)

        # XXX keep as a function (we may need hexagonal grids later!)
        X, Y, Z = self.cube_meshgrid(
            x_lims, y_lims, z_lims, nx, ny, nz
        )

        surface = f(X, Y, Z, **kwargs)
        vertices, faces = self._extract_features(surface, level)

        self._vertices = vertices
        self._faces    = faces
        self._mesh     = Trimesh(vertices=vertices, faces=faces)

        # 1. Map the raw surface mesh to physical space first!
        self._mesh = self.map_voxel_to_physical(
            mesh   = self._mesh,
            x_lims = x_lims,
            y_lims = y_lims,
            z_lims = z_lims,
            nx     = nx,
            ny     = ny,
            nz     = nz,
        )

        # 2. Thicken the physical surface mesh using physical thickness!
        if thickness > 0.0:
            thickened = self.thicken_surface_mesh(
                self._mesh, thickness
            )
            self._mesh = thickened

        # 3. Always fix normals and volume if it is a closed shell!
        self._mesh.fix_normals()

        if self._mesh.is_volume and self._mesh.volume < 0:
            self._mesh.invert()

        self._vertices = self._mesh.vertices
        self._faces    = self._mesh.faces

    def _resolve_functional(
            self,
            functional: str | FunctionalType,
        ) -> FunctionalType:
        """ Resolve a functional to be used for mesh generation. """
        if callable(functional):
            return functional

        if isinstance(functional, str):
            if not (module := Path(functional)).exists():
                raise FileNotFoundError(f"No such module {functional}")

            return self.load_implicit_surface(module)

        raise ValueError("Functional must be a string or a callable")

    def _extract_features(
            self,
            surface: NDArray[Any],
            fractional_level: float,
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
    def load_implicit_surface(path: Path) -> FunctionalType:
        """ Import ``implicit_surface``  from path to module.

        Parameters
        ----------
        path : Path
            Path to a Python file defining ``implicit_surface``.

        Raises
        ------
        FileNotFoundError
            If path does not point to an existing file.
        AttributeError
            If the loaded module does not expose ``implicit_surface``.
        """
        if not (resolved := path.resolve()).is_file():
            raise FileNotFoundError(f"No such file: {resolved}")

        spec = importlib.util.spec_from_file_location(
            resolved.stem, resolved
        )

        if spec is None or spec.loader is None:
            raise ImportError(f"Cannot load module from: {resolved}")

        module: ModuleType = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)  # type: ignore[union-attr]

        if not hasattr(module, "implicit_surface"):
            raise AttributeError(
                f"Module '{resolved}' does not define 'implicit_surface'"
            )

        return module.implicit_surface  # type: ignore[return-value]

    @staticmethod
    def cube_meshgrid(
            x_lims: tuple[float, float],
            y_lims: tuple[float, float],
            z_lims: tuple[float, float],
            nx: int,
            ny: int,
            nz: int,
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

    @staticmethod
    def cube_subvolume(
            stl_mesh: Trimesh,
            box_frac: float = 0.98,
        ) -> Trimesh:
        """ Create a region that is a fraction of the bounding box.

        Parameters
        ----------
        stl_mesh : trimesh.Trimesh
            The STL mesh of the gyroid.
        box_frac : float, default=0.98
            Fraction of the bounding box dimensions to scale by.

        Returns
        -------
        trimesh.Trimesh
            Sub-volume box centered on the original bounding box.
        """
        x_lims, y_lims, z_lims = stl_mesh.bounding_box.bounds.T

        cx = (x_lims[0] + x_lims[1]) / 2.0
        cy = (y_lims[0] + y_lims[1]) / 2.0
        cz = (z_lims[0] + z_lims[1]) / 2.0

        dx = (x_lims[1] - x_lims[0]) * box_frac
        dy = (y_lims[1] - y_lims[0]) * box_frac
        dz = (z_lims[1] - z_lims[0]) * box_frac

        box = trimesh.creation.box(extents=(dx, dy, dz))
        box.apply_translation([cx, cy, cz])

        return box

    @staticmethod
    def thicken_surface_mesh(
            mesh: Trimesh,
            thickness: float,
            *,
            centered: bool = True,
        ) -> Trimesh:
        """ Create a shell with finite thickness from a surface.

        The algorithm offsets vertices along vertex normals to form
        two surfaces and connects boundary edges with side triangles.

        Parameters
        ----------
        mesh : Trimesh
            Input surface mesh.
        thickness : float
            Total shell thickness (must be positive).
        centered : bool, default=True
            If true, offsets equally to both sides of the surface.
            If false, keeps original surface as outer side.

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
        faces    = np.asarray(mesh.faces, dtype=np.int32)
        normals  = np.asarray(mesh.vertex_normals)

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

        # Find all directed boundary edges to enforce winding
        edges = np.vstack((
            faces[:, [0, 1]],
            faces[:, [1, 2]],
            faces[:, [2, 0]],
        ))

        sorted_edges = np.sort(edges, axis=1)

        unique_edges, counts = np.unique(
            sorted_edges, axis=0, return_counts=True
        )

        boundary_undirected = set(
            tuple(e) for e in unique_edges[counts == 1]
        )

        boundary_directed = []

        for face in faces:
            pairs = [
                (face[0], face[1]),
                (face[1], face[2]),
                (face[2], face[0]),
            ]

            for u, v in pairs:
                if tuple(sorted((u, v))) in boundary_undirected:
                    boundary_directed.append((u, v))

        if len(boundary_directed) > 0:
            boundary_directed = np.array(boundary_directed, dtype=np.int32)

            u = boundary_directed[:, 0]
            v = boundary_directed[:, 1]

            u_in = u + n_vertices
            v_in = v + n_vertices

            side_1 = np.column_stack((v, u, u_in))
            side_2 = np.column_stack((v, u_in, v_in))

            side_faces = np.vstack((side_1, side_2)).astype(np.int32)
            shell_faces = np.vstack((outer_faces, inner_faces, side_faces))
        else:
            shell_faces = np.vstack((outer_faces, inner_faces))

        shell_vertices = np.vstack((outer_vertices, inner_vertices))
        return Trimesh(
            vertices=shell_vertices, faces=shell_faces, process=True
        )

    @staticmethod
    def map_voxel_to_physical(
            mesh: Trimesh,
            x_lims: tuple[float, float],
            y_lims: tuple[float, float],
            z_lims: tuple[float, float],
            nx: int,
            ny: int,
            nz: int,
        ) -> Trimesh:
        """ Scale mesh vertices from voxel to physical coordinates.

        Marching cubes extracts vertices in index space [0, N-1],
        which needs to be mapped back to the physical ranges.

        Note: We account for the 'xy' indexing of np.meshgrid in
        FunctionalShapes.cube_meshgrid, which swaps the X and Y axes.
        """
        pts = mesh.vertices.copy()
        pts_physical = np.zeros_like(pts)

        delta_x = (x_lims[1] - x_lims[0]) / (nx - 1)
        delta_y = (y_lims[1] - y_lims[0]) / (ny - 1)
        delta_z = (z_lims[1] - z_lims[0]) / (nz - 1)

        pts_physical[:, 0] = x_lims[0] + pts[:, 1] * delta_x
        pts_physical[:, 1] = y_lims[0] + pts[:, 0] * delta_y
        pts_physical[:, 2] = z_lims[0] + pts[:, 2] * delta_z

        mesh.vertices = pts_physical

        # Fix chiral reflection: a coordinate swap inverts the winding of
        # all faces. We must fix normals and invert faces if volume is
        # negative to keep the mesh as a solid volume.
        mesh.fix_normals()

        if mesh.is_volume and mesh.volume < 0:
            mesh.invert()

        return mesh

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

    @staticmethod
    def from_stl(filename: str | Path, **kwargs: Any) -> Trimesh:
        """ Load an STL mesh file.

        Note: if the STL was not created via FunctionalShapes, you may
        need to call `FunctionalShapes.map_voxel_to_physical` on the
        mesh from STL to get the right voxel mapping.
        """
        mesh = trimesh.load(filename, file_type="stl")

        if kwargs.pop("map_voxel_to_physical", False):
            mesh = FunctionalShapes.map_voxel_to_physical(mesh, **kwargs)

        return mesh

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


class PorousDomainExtractor:
    """ Extract fluid and solid domains from STL mesh and bounding volume.

    This class computes the fluid (pore) domains by subtracting the
    thickened surface from the bounding sub-volume, and computes
    the solid domains by calculating the intersection between them.

    Parameters
    ----------
    mesh : trimesh.Trimesh
        The thickened physical solid surface mesh.
    volume : trimesh.Trimesh
        The bounding sub-volume mesh.
    """

    __slots__ = (
        "_mesh",
        "_volume",
        "_bodies",
        "_fluid_meshes",
        "_solid_meshes",
    )

    def __init__(self,
            mesh: trimesh.Trimesh,
            volume: trimesh.Trimesh
        ) -> None:
        self._mesh   = mesh
        self._volume = volume
        self._bodies = self.get_porous_domain(mesh, volume)

        self._fluid_meshes: list[trimesh.Trimesh] = []
        self._solid_meshes: list[trimesh.Trimesh] = []

    @staticmethod
    def _filter_bodies(
            bodies: list[trimesh.Trimesh],
            min_vertices: int = 100,
        ) -> list[trimesh.Trimesh]:
        """ Filter out small bodies and fix their orientation.

        Parameters
        ----------
        bodies : list of trimesh.Trimesh
            The list of bodies to filter.
        min_vertices : int, default=100
            Minimum number of vertices for a body to be retained.

        Returns
        -------
        list of trimesh.Trimesh
            The filtered and corrected list of bodies.
        """
        result = []

        for c in bodies:

            if len(c.vertices) < min_vertices:
                continue

            c.fix_normals()

            if c.volume < 0:
                c.invert()

            result.append(c)

        return result

    @staticmethod
    def export_stl_ascii(
            mesh: trimesh.Trimesh,
            filepath: str | Path,
        ) -> None:
        """ Export a trimesh to an ASCII STL file.

        The first line of the solid name matches the filename.

        Parameters
        ----------
        mesh : trimesh.Trimesh
            The trimesh object to export.
        filepath : str or Path
            The file path where the STL file will be saved.
        """
        name = Path(filepath).stem
        ascii_stl_string = trimesh.exchange.stl.export_stl_ascii(mesh)
        lines = ascii_stl_string.splitlines()

        if lines:

            if lines[0].startswith("solid"):
                lines[0] = f"solid {name}"

            if lines[-1].startswith("endsolid"):
                lines[-1] = f"endsolid {name}"

        ascii_stl_string = "\n".join(lines) + "\n"

        with open(filepath, "w", encoding="utf-8") as f:
            f.write(ascii_stl_string)

    @classmethod
    def _operation(
            cls,
            this: trimesh.Trimesh,
            that: trimesh.Trimesh,
            operation: str,
        ) -> trimesh.Trimesh:
        """ Perform robust CSG boolean operations via Blender.

        Parameters
        ----------
        this : trimesh.Trimesh
            First operand mesh.
        that : trimesh.Trimesh
            Second operand mesh.
        operation : str
            The name of the trimesh boolean operation.

        Returns
        -------
        trimesh.Trimesh
            The resulting boolean mesh.
        """
        options = dict(engine="blender", check_volume=False, use_exact=False)
        return getattr(this, operation)(that, **options)

    @property
    def bodies(self) -> list[trimesh.Trimesh]:
        """ Return all extracted fluid and solid bodies. """
        return self._bodies

    @property
    def solid_bodies(self) -> list[trimesh.Trimesh]:
        """ Return all extracted solid bodies. """
        if not self._solid_meshes:
            self._solid_meshes = [
                b for b in self.bodies if b.metadata["type"] == "solid"
            ]

        return self._solid_meshes

    @property
    def fluid_bodies(self) -> list[trimesh.Trimesh]:
        """ Return all extracted fluid bodies. """
        if not self._fluid_meshes:
            self._fluid_meshes = [
                b for b in self.bodies if b.metadata["type"] == "fluid"
            ]

        return self._fluid_meshes

    @classmethod
    def get_porous_domain(
            cls,
            stl_physical: trimesh.Trimesh,
            stl_domain: trimesh.Trimesh,
            min_vertices: int = 100,
        ) -> list[trimesh.Trimesh]:
        """ Extract and split fluid and solid components.

        Parameters
        ----------
        stl_physical : trimesh.Trimesh
            The physical solid wall mesh.
        stl_domain : trimesh.Trimesh
            The bounding domain box mesh.
        min_vertices : int, default=100
            Minimum vertex threshold for filtering.

        Returns
        -------
        list of trimesh.Trimesh
            The partitioned fluid and solid components.
        """
        # Determine appropriate scale factor to avoid Blender tolerance issues
        extents = stl_domain.extents
        max_extent = float(np.max(extents))
        scale_factor = 100.0 / max_extent if max_extent > 0.0 else 1.0

        if scale_factor != 1.0:
            stl_domain_scaled = stl_domain.copy()
            stl_domain_scaled.apply_scale(scale_factor)
            stl_physical_scaled = stl_physical.copy()
            stl_physical_scaled.apply_scale(scale_factor)
        else:
            stl_domain_scaled = stl_domain
            stl_physical_scaled = stl_physical

        def by_vertices(c):
            return len(c.vertices)

        # 1. Fluid domain: subtract the solid wall from the sub-volume
        pore_trimesh = cls._operation(
            stl_domain_scaled, stl_physical_scaled, "difference"
        )
        porous = pore_trimesh.split(only_watertight=False)
        porous = cls._filter_bodies(porous, min_vertices=min_vertices)

        for c in porous:
            c.metadata["type"] = "fluid"

        # 2. Solid domain: intersect the solid wall with the sub-volume
        solid_trimesh = cls._operation(
            stl_domain_scaled, stl_physical_scaled, "intersection"
        )
        solids = solid_trimesh.split(only_watertight=False)
        solids = cls._filter_bodies(solids, min_vertices=min_vertices)

        for c in solids:
            c.metadata["type"] = "solid"

        # Scale all bodies back down to original coordinate space
        if scale_factor != 1.0:
            inv_scale = 1.0 / scale_factor

            for c in porous:
                c.apply_scale(inv_scale)

            for c in solids:
                c.apply_scale(inv_scale)

        # Sort largest components first
        porous = sorted(porous, key=by_vertices, reverse=True)
        solids = sorted(solids, key=by_vertices, reverse=True)

        return porous + solids

    def save_project(
            self,
            project_dir: Path,
            overwrite: bool = False,
            save_sources: bool = False,
        ) -> None:
        """ Save extracted domains to a target directory.

        Parameters
        ----------
        project_dir : Path
            The directory path where STL outputs will be written.
        overwrite : bool, default=False
            Whether to overwrite the directory if it exists.
        save_sources : bool, default=False
            Whether to also save the original input source meshes.

        Raises
        ------
        FileExistsError
            If the project directory exists and overwrite is False.
        """

        if project_dir.exists() and not overwrite:
            raise FileExistsError(f"{project_dir} already exists.")

        project_dir.mkdir(parents=True, exist_ok=True)

        if save_sources:
            self.export_stl_ascii(
                self._mesh, project_dir / "source_mesh.stl"
            )
            self.export_stl_ascii(
                self._volume, project_dir / "source_volume.stl"
            )

        for i, body in enumerate(self.fluid_bodies):
            self.export_stl_ascii(body, project_dir / f"fluid_{i}.stl")

        for i, body in enumerate(self.solid_bodies):
            self.export_stl_ascii(body, project_dir / f"solid_{i}.stl")


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
    """ Decorator to apply noise and sinusoidal modulation. """
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

        surface = noisy_surface(surface, A * abs(vmax - vmin))
        vmin, vmax = surface.min(), surface.max()

        surface = wavy_surface(surface, Z, B * abs(vmax - vmin), f)

        return surface

    return decorated


def plot_domain(
        bodies: list[Trimesh] | None = None,
        volume: Trimesh | None = None,
        parent: Trimesh | None = None,
        saveas: str | Path | None = None,
    ) -> None:
    """ Display the extracted domain components.

    Parameters
    ----------
    bodies : list of trimesh.Trimesh, optional
        List of extracted domain components (fluid and solid).
    volume : trimesh.Trimesh, optional
        Sub-volume box mesh.
    parent : trimesh.Trimesh, optional
        Original solid domain mesh.
    saveas : str or Path, optional
        File path to save the generated rendering as a PNG.
    """
    off_screen = saveas is not None
    plotter = pv.Plotter(off_screen=off_screen, window_size=[800, 600])

    # Add thin reference model of the solid wall shell
    if parent is not None:
        solid_pv = pv.wrap(parent)
        plotter.add_mesh(
            mesh        = solid_pv,
            color       = "#888888",
            opacity     = 1.0,
            show_edges  = False,
            label       = "Solid Wall (Reference)"
        )

    # Add wireframe bounding box.
    if volume is not None:
        solid_pv = pv.wrap(volume)
        plotter.add_mesh(
            solid_pv.outline(), color="#000000", line_width=2
        )

    has_labels = False

    if bodies is not None:
        fluid_colors = ["#012169", "#FF8200", "#A50033", "#009681"]
        fluid_idx = 0
        solid_idx = 0

        for body in bodies:
            body_type = body.metadata.get("type", "fluid")

            if body_type == "fluid":
                color = fluid_colors[fluid_idx % len(fluid_colors)]
                label = (
                    f"Fluid Zone {fluid_idx + 1} "
                    f"(Vol: {body.volume:.4f})"
                )
                opacity = 0.5
                fluid_idx += 1
            else:
                color = "#7F7F7F"
                label = (
                    f"Solid Zone {solid_idx + 1} "
                    f"(Vol: {body.volume:.4f})"
                )
                opacity = 1.0
                solid_idx += 1

            plotter.add_mesh(
                mesh        = body,
                color       = color,
                opacity     = opacity,
                show_edges  = False,
                label       = label
            )
            has_labels = True

    if has_labels or parent is not None:
        plotter.add_legend(face="circle", bcolor=None)

    plotter.show_grid(
        xtitle="X [m]",
        ytitle="Y [m]",
        ztitle="Z [m]",
        color="#555555",
        font_size=10,
        location="outer",
        fmt="{:.2e}",
    )
    plotter.show_axes()

    if off_screen:
        plotter.screenshot(saveas)
    else:
        plotter.show()


def main() -> int:
    #region: build parser
    parser = argparse.ArgumentParser(
        prog="porous-media-geometry",
        description="Generate porous media geometry from YAML config.",
    )

    parser.add_argument(
        "config",
        type=Path,
        help="Path to YAML file with FunctionalShapes parameters.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting existing output files.",
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

# Geometry

```python
from skimage import measure
from trimesh import Trimesh
import matplotlib.pyplot as plt
import numpy as np
```

```python
class MeshGrid3D:
    def __init__(self, x_rng, y_rng, z_rng, nx, ny, nz):
        self._x = np.linspace(*x_rng, nx)
        self._y = np.linspace(*y_rng, ny)
        self._z = np.linspace(*z_rng, nz)

    def get_scaled(self, sx=1, sy=1, sz=1):
        x = self._x * sx
        y = self._y * sy
        z = self._z * sz
        scale = np.array([sx, sy, sz])
        grid = np.meshgrid(x, y, z)
        return grid, scale


class Shapes:
    @staticmethod
    def gyroid(X, Y, Z):
        a = np.sin(X) * np.cos(Y)
        b = np.sin(Y) * np.cos(Z)
        c = np.sin(Z) * np.cos(X)
        return a + b + c

    @staticmethod
    def schwarz(X, Y, Z, isocontour=0.0):
        return np.cos(X) + np.cos(Y) + np.cos(Z) - isocontour

    @staticmethod
    def to_surface(f, rugosity=0.0, stripes=0.0, freq=1.0):
        rugosity = rugosity * np.random.normal(size=f.shape)
        stripes  = stripes * np.sin(freq * Z)
        return f + stripes + rugosity

    @staticmethod
    def extract_features(surface, level=0):
        the_min, the_max = extrema(surface)

        if not (the_min <= level <= the_max):
            raise ValueError(f"Level {level} outside [{the_min}; {the_max}]")
            
        verts, faces, _, _ = measure.marching_cubes(surface, level=level)
        faces = faces.astype(np.int32)

        if len(verts) == 0 or len(faces) == 0:
            raise RuntimeError("No surface generated!")
        
        return verts, faces


def plot_trimesh(verts, faces, scale, **kwargs):
    vx = verts[:, 0] * scale[0]
    vy = verts[:, 1] * scale[1]
    vz = verts[:, 2] * scale[2]

    opts = dict(
        cmap  = kwargs.get('cmap', 'gnuplot'), 
        lw    = kwargs.get('lw', 0.1), 
        alpha = kwargs.get('alpha', 1.0)
    )
    
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='3d')
    ax.plot_trisurf(vx, vy, vz, triangles=faces, **opts)
    ax.auto_scale_xyz(extrema(vx), extrema(vy), extrema(vz))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    fig.tight_layout()

    return fig, ax


def extrema(a):
    return [a.min(), a.max()]
```

```python
grid = MeshGrid3D(
    x_rng = (-3*np.pi, 6*np.pi),
    y_rng = (-3*np.pi, 6*np.pi),
    z_rng = (       0, 3*np.pi),
    nx    = 50,
    ny    = 50,
    nz    = 50
)
(X, Y, Z), scale = grid.get_scaled(sx=1, sy=1, sz=1)

# surf = Shapes.gyroid(X, Y, Z, isocontour=0.1)
surface = Shapes.schwarz(X, Y, Z, isocontour=0.1)
surface = Shapes.to_surface(surface, rugosity=0.0, stripes=0.0, freq=5.0)
verts, faces = Shapes.extract_features(surface, level=0.5)

mesh = Trimesh(vertices=verts, faces=faces)
_ = mesh.export(file_obj='surface.stl')

_ = plot_trimesh(verts, faces, scale)
```

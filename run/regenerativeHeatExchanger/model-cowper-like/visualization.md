# Results visualization

```python
import pyvista as pv
from heat_recovery.cowper_like import CowperLikePost
from heat_recovery.thermocline import ThermoclineModel

pv.set_jupyter_backend("static")
```

```python
model = ThermoclineModel("dimensioning.yaml")
```

```python
post = CowperLikePost(scale=(1, 1, model.num_D_h / model.num_h_t))
post.load_state()
```

```python
opts = post.get_options()
pl = pv.Plotter()

pl.add_mesh(post.slice_fluid, scalars="T", cmap="hot", **opts)
pl.add_mesh(post.slice_solid, scalars="T", cmap="hot", show_scalar_bar=False)
CowperLikePost.align_camera(pl, xc=0.025, zc=0.02, ps=0.017)

pl.show()
```

```python
opts = post.get_options()
pl = pv.Plotter(shape=(1, 2), border=False)

pl.subplot(0, 0)
pl.add_mesh(post.slice_fluid, scalars="pRel", cmap="jet", **opts)
CowperLikePost.align_camera(pl, xc=0.0125, zc=0.02, ps=0.017)

pl.subplot(0, 1)
pl.add_mesh(post.slice_fluid, scalars='p_rgh', cmap="jet", **opts)
CowperLikePost.align_camera(pl, xc=0.0125, zc=0.02, ps=0.017)

pl.show()
```

```python
opts = post.get_options()
pl = pv.Plotter(shape=(1, 2), border=False)

pl.subplot(0, 0)
pl.add_mesh(post.slice_fluid, scalars='pRel', cmap="jet", **opts)
CowperLikePost.align_camera(pl, xc=0.0125, zc=0.02, ps=0.017)

pl.subplot(0, 1)
pl.add_mesh(post.slice_fluid, scalars='rho', cmap="jet", **opts)
CowperLikePost.align_camera(pl, xc=0.0125, zc=0.02, ps=0.017)

pl.show()
```

```python
opts = post.get_options()
pl = pv.Plotter()

# pl.add_mesh(post.slice_fluid, scalars="Umag", cmap="jet", **opts)
pl.add_mesh(post.slice_fluid, scalars="U", component=2, cmap="jet", **opts)
CowperLikePost.align_camera(pl, xc=0.0125, zc=0.02, ps=0.017)

pl.show()
```

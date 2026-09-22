#import "../book.typ": *

= Inflation layers <chapter-4>

In computational fluid dynamics, accurately predicting aerodynamic forces, wall shear stress, skin friction, heat transfer, and flow separation requires capturing the steep velocity and thermal gradients present in boundary layers. Standard Cartesian volume cells generated during castellation and snapping, while mathematically valid, are often too isotropic or too coarse near boundaries to resolve these thin shear layers without requiring an excessive global cell count. Boundary layer inflation introduces thin, anisotropic, prismatic cell layers aligned with solid boundary surfaces. While volume mesh resolution in `snappyHexMesh` is constrained by a strict two-to-one octree division ratio, the layer generation routine provides precise control over individual layer heights, expansion factors, and total layer counts along specific boundaries.

== The layer addition process

The inflation layer generation procedure is governed by the `addLayersControls` sub-dictionary within `snappyHexMeshDict` and activated via the `addLayers` top-level boolean switch. Unlike volume cell splitting, layer addition is an extrusion and morphing process executed directly upon the body-fitted snapped mesh. The algorithm operates through an iterative sequence:

- The boundary vertices of the snapped mesh are displaced inward, projecting backward into the fluid domain along surface normal vectors by a designated total layer thickness.

- A displacement relaxation equation is solved across the interior mesh to accommodate this inward boundary displacement without corrupting internal cell topologies.

- Geometric mesh quality metrics are evaluated on the distorted internal cells; if quality criteria are violated, the projected extrusion distance is reduced iteratively.

- Once validation criteria are satisfied across the domain, prismatic cell layers are inserted into the created void between the retracted volume mesh and the surface boundary.

- A comprehensive post-insertion quality check is executed; if any inserted or neighboring cells fail mesh quality constraints, the conflicting local layers are removed or collapsed back to maintain overall numerical validity.

== Layer sizing strategies

Within `addLayersControls`, the geometric distribution of the inflation stack is specified using the `layers` sub-dictionary and a set of layer thickness keywords. The `layers` sub-dictionary explicitly declares the target surface patches and their required integer layer counts using the `nSurfaceLayers` parameter. Because layer addition operates on the active mesh rather than the source CAD geometry, these entries reference the patch names generated during the snapping stage rather than raw triangulated surface regions.

To define the dimensional progression of the prismatic cells, `snappyHexMesh` provides four primary thickness parameters: `expansionRatio`, `firstLayerThickness`, `finalLayerThickness`, and total `thickness`. The user must specify exactly two of these four parameters. Supplying more than two parameters creates an over-constrained mathematical system, while fewer results in an under-determined setup. The `relativeSizes` boolean switch dictates whether these dimensions are interpreted as absolute physical distances or as fractional ratios relative to the local undistorted volume cell size directly adjacent to the boundary layer. When configuring wall-resolved turbulence simulations, combining an absolute `firstLayerThickness` with an `expansionRatio` allows precise targeting of non-dimensional wall distance values, whereas combining `finalLayerThickness` and `expansionRatio` with relative sizing ensures a smooth geometric transition between the outermost prism layer and the adjacent background hexahedral cells.

== Topological controls and surface smoothing

// TODO create an illustration of the featureAngle definition/process
Near-wall mesh quality during extrusion is highly sensitive to surface curvature, sharp corners, and intersecting features. To prevent distorted, concave, or self-intersecting prismatic cells, `addLayersControls` includes several topological and geometric control parameters. The `featureAngle` setting specifies the maximum geometric angle across which layers are allowed to continuously extrude; along sharp edges exceeding this threshold, layer growth is terminated to prevent overlapping normal vectors. The `nGrow` parameter governs the number of connected cell faces adjacent to non-extruded points that are progressively stepped down, avoiding abrupt cliff-like terminations near complex topological features.

To maintain uniform layer thickness across curved boundaries, surface normals and interior movement vectors are smoothed iteratively using `nSmoothSurfaceNormals` and `nSmoothNormals`, while `nSmoothThickness` averages overall layer thickness across adjacent surface faces. Geometric distortion is further regulated by `maxFaceThicknessRatio`, which halts layer extrusion across heavily warped faces, and `maxThicknessToMedialRatio`, which reduces layer thickness in tight internal corners or narrow passages where opposing boundary layers approach one another along the medial axis.

== Troubleshooting and iterative strategies

Layer addition is notoriously recognized as the most challenging and sensitive phase of the `snappyHexMesh` pipeline. In complex geometries with sharp trailing edges or narrow gaps, default parameters frequently lead to partial layer collapse or localized deletion due to quality check failures.

To diagnose and resolve these issues efficiently, users should adopt a modular workflow. By keeping `addLayers` set to false while finalizing `castellatedMesh` and `snap`, the core body-fitted mesh can be fully validated before initiating layer inflation. Once a high-quality snapped mesh is achieved, the case can be restarted from the snapped time directory with `castellatedMesh` and `snap` deactivated and `addLayers` set to true. When layer insertion fails to reach the requested layer count, relaxing internal quality constraints via the `relaxed` sub-dictionary within `meshQualityControls`, increasing `nLayerIter` and `nRelaxIter`, or locally increasing surface refinement to yield smaller, more flexible base cells will significantly improve layer coverage and mesh validity.

== Add layers controls parameters

=== Layer sizing and dimensions

Layer sizing controls determine prism height and expansion outward from the boundary surface. The interpretation of these values is governed by `relativeSizes`, which is a boolean switch (`true` or `false`). If active, thicknesses are specified as fractions of the local undistorted background cell height directly adjacent to the patch. Setting it to `false`, thicknesses are defined as absolute physical dimensions in domain units (typically meters).

#exampleblock(title: "the relative size of a cell")[
If setting `relativeSizes true`, then if an adjacent background volume cell has an edge height $Delta x = 2.0 "mm"$, setting a thickness of `0.5` equates to an absolute thickness of $1.0 "mm"$.
]

#exampleblock(title: "the absolute size of a cell")[
  If setting `relativeSizes false`, then `firstLayerThickness 0.0005` creates layers with physical dimensions starting at $0.5 "mm"$, regardless of background mesh refinement.
]

Users must specify *exactly two* of the following four geometric progression parameters (the remaining two are calculated automatically to avoid over-constraining the geometric series):

- `expansionRatio`: The geometric expansion factor between consecutive layers moving away from the wall. A value of `1.2` means each successive layer outward is 20% thicker than the preceding layer.

- `finalLayerThickness`: The target thickness of the outermost layer (the layer adjacent to the internal volume mesh). When `relativeSizes` is `true`, this sets the outermost layer to a given fraction of the adjacent undistorted background cell height to preserve a smooth volumetric size transition.

- `firstLayerThickness`:

- `thickness`:

Other layer sizing parameters:

- `minThickness`: The minimum permissible thickness for an individual cell layer (relative to the background cell or absolute, depending on `relativeSizes`). If mesh shrinking or quality-optimization steps compress a layer below this value, layer addition is canceled or collapsed locally.

=== Extrusion topology and features

- `featureAngle`: The maximum allowable angle between adjacent surface face normals across which layers can continuously extrude. In OpenFOAM, 0° corresponds to a flat plane and 90° represents perpendicular faces. If the normal turning angle across an edge exceeds 100°, layer extrusion terminates across that edge to prevent self-intersecting displacement vectors.

- `slipFeatureAngle`: Controls sliding at intersecting, non-extruded boundary patches (such as symmetry planes, inlets, or slip walls). If the angle between the layer extrusion direction and the intersecting patch normal is greater than `slipFeatureAngle` (typically defaulting to $0.5 times "featureAngle"$), the mesh displacement is permitted to slip along that patch instead of sticking or collapsing rigidly.

- `nGrow`: Specifies how many additional rings of connected faces around non-extruded points are flagged to terminate layer growth. Setting `nGrow 0` allows layers to terminate immediately at the offending feature. A positive integer grows the non-extrusion zone wider, helping mesh relaxation and convergence near sharp corners by providing a buffer against pinching.

- `nBufferCellsNoExtrude`: Specifies the number of buffer cells used to step down layer thickness and count towards terminating edges. Rather than ending a 5-layer stack abruptly against a wall or non-extruded face, positive values create a gradual step-down transition (_e.g._, tapering from 5 to 4, 3, 2, 1 layers) to avoid abrupt cell mismatches in the volume mesh.

=== Iteration and Convergence Controls

- `nLayerIter`: The maximum overall iterations allowed for the layer addition process. snappyHexMesh will attempt up to this value in number of cycles of projecting the boundary vertices outward, smoothing internal displacements, and validating mesh quality.

- `nRelaxedIter`: The iteration threshold beyond which the mesh mover switches from the strict quality criteria to the relaxed criteria defined in the `relaxed` sub-dictionary of `meshQualityControls`. If layer insertion struggles to converge within 20 iterations under standard metrics (e.g., non-orthogonality limits), relaxed settings are applied to encourage completion.

=== Smoothing controls

- `nSmoothSurfaceNormals`: The number of smoothing sweeps applied directly to the surface normal vectors along the boundary patches before computing initial layer extrusion paths.

- `nSmoothNormals`: The number of smoothing sweeps applied to the interior mesh movement direction vectors. This smears the displacement directions into the internal domain, mitigating face skewness and edge collisions as the internal mesh compresses.

- `nSmoothThickness`: The number of smoothing iterations applied to the layer thickness field across patch faces. This produces a continuous, gradual variation in layer thickness across adjoining faces and patches.

=== Medial axis and quality limiting

- `maxFaceThicknessRatio`: Terminates layer growth on boundary faces where the total extrusion thickness exceeds 50% of the face's characteristic size (or where faces are excessively warped). This prevents inverted or high-aspect-ratio degenerate cells along severely stretched boundary elements.

- `maxThicknessToMedialRatio`: Truncates layer thickness in narrow passages, acute corners, or internal gaps where opposing walls converge. snappyHexMesh identifies the topological centerline between opposing surfaces (the medial axis). If total layer thickness exceeds this fraction of the distance to the medial axis, layer thickness is scaled down to avoid self-collision and negative cell volumes.

- `minMedianAxisAngle`: The angular criterion used to detect the medial axis / centerline. When the vectors pointing toward nearest wall features differ across an edge by an angle sharper than 90°, that edge is designated as a medial axis location for the distance-field wave solver.

- `nMedialAxisIter`: The maximum number of sweeps or diffusion iterations allocated to propagate, calculate, and smooth the medial axis distance field across the mesh domain.

// EOF
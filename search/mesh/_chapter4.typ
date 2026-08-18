= Inflation layers <chapter-4>

// Physical Significance and Near-Wall Resolution

In computational fluid dynamics, accurately predicting aerodynamic forces, wall shear stress, skin friction, heat transfer, and flow separation requires capturing the steep velocity and thermal gradients present in boundary layers. Standard Cartesian volume cells generated during castellation and snapping, while mathematically valid, are often too isotropic or too coarse near boundaries to resolve these thin shear layers without requiring an excessive global cell count. Boundary layer inflation introduces thin, anisotropic, prismatic cell layers aligned with solid boundary surfaces. While volume mesh resolution in `snappyHexMesh` is constrained by a strict two-to-one octree division ratio, the layer generation routine provides precise control over individual layer heights, expansion factors, and total layer counts along specific boundaries.

// Mechanics of the Layer Addition Routine

The inflation layer generation procedure is governed by the `addLayersControls` sub-dictionary within `snappyHexMeshDict` and activated via the `addLayers` top-level boolean switch. Unlike volume cell splitting, layer addition is an extrusion and morphing process executed directly upon the body-fitted snapped mesh. The algorithm operates through an iterative sequence:

The boundary vertices of the snapped mesh are displaced inward, projecting backward into the fluid domain along surface normal vectors by a designated total layer thickness.

A displacement relaxation equation is solved across the interior mesh to accommodate this inward boundary displacement without corrupting internal cell topologies.

Geometric mesh quality metrics are evaluated on the distorted internal cells; if quality criteria are violated, the projected extrusion distance is reduced iteratively.

Once validation criteria are satisfied across the domain, prismatic cell layers are inserted into the created void between the retracted volume mesh and the surface boundary.

A comprehensive post-insertion quality check is executed; if any inserted or neighboring cells fail mesh quality constraints, the conflicting local layers are removed or collapsed back to maintain overall numerical validity.

// Layer Sizing Strategies and Mathematical Formulation

Within `addLayersControls`, the geometric distribution of the inflation stack is specified using the `layers` sub-dictionary and a set of layer thickness keywords. The `layers` sub-dictionary explicitly declares the target surface patches and their required integer layer counts using the `nSurfaceLayers` parameter. Because layer addition operates on the active mesh rather than the source CAD geometry, these entries reference the patch names generated during the snapping stage rather than raw triangulated surface regions.

To define the dimensional progression of the prismatic cells, `snappyHexMesh` provides four primary thickness parameters: `expansionRatio`, `firstLayerThickness`, `finalLayerThickness`, and total `thickness`. The user must specify exactly two of these four parameters. Supplying more than two parameters creates an over-constrained mathematical system, while fewer results in an under-determined setup. The `relativeSizes` boolean switch dictates whether these dimensions are interpreted as absolute physical distances or as fractional ratios relative to the local undistorted volume cell size directly adjacent to the boundary layer. When configuring wall-resolved turbulence simulations, combining an absolute `firstLayerThickness` with an `expansionRatio` allows precise targeting of non-dimensional wall distance values, whereas combining `finalLayerThickness` and `expansionRatio` with relative sizing ensures a smooth geometric transition between the outermost prism layer and the adjacent background hexahedral cells.

// Topological Controls and Surface Smoothing Parameters

Near-wall mesh quality during extrusion is highly sensitive to surface curvature, sharp corners, and intersecting features. To prevent distorted, concave, or self-intersecting prismatic cells, `addLayersControls` includes several topological and geometric control parameters. The `featureAngle` setting specifies the maximum geometric angle across which layers are allowed to continuously extrude; along sharp edges exceeding this threshold, layer growth is terminated to prevent overlapping normal vectors. The `nGrow` parameter governs the number of connected cell faces adjacent to non-extruded points that are progressively stepped down, avoiding abrupt cliff-like terminations near complex topological features.

To maintain uniform layer thickness across curved boundaries, surface normals and interior movement vectors are smoothed iteratively using `nSmoothSurfaceNormals` and `nSmoothNormals`, while `nSmoothThickness` averages overall layer thickness across adjacent surface faces. Geometric distortion is further regulated by `maxFaceThicknessRatio`, which halts layer extrusion across heavily warped faces, and `maxThicknessToMedialRatio`, which reduces layer thickness in tight internal corners or narrow passages where opposing boundary layers approach one another along the medial axis.

// Troubleshooting and Iterative Workflow Strategies

Layer addition is notoriously recognized as the most challenging and sensitive phase of the `snappyHexMesh` pipeline. In complex geometries with sharp trailing edges or narrow gaps, default parameters frequently lead to partial layer collapse or localized deletion due to quality check failures.

To diagnose and resolve these issues efficiently, users should adopt a modular workflow. By keeping `addLayers` set to false while finalizing `castellatedMesh` and `snap`, the core body-fitted mesh can be fully validated before initiating layer inflation. Once a high-quality snapped mesh is achieved, the case can be restarted from the snapped time directory with `castellatedMesh` and `snap` deactivated and `addLayers` set to true. When layer insertion fails to reach the requested layer count, relaxing internal quality constraints via the `relaxed` sub-dictionary within `meshQualityControls`, increasing `nLayerIter` and `nRelaxIter`, or locally increasing surface refinement to yield smaller, more flexible base cells will significantly improve layer coverage and mesh validity.

// EOF
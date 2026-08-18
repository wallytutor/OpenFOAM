= Mesh quality and validation <chapter-5>

// Paradigms of Spatial Grid Evaluation

In computational fluid dynamics and multiphysics simulations, the spatial grid establishes the discrete mathematical framework required to solve the governing conservation equations. Achieving a high-quality mesh is an essential prerequisite for obtaining accurate, stable, and physically realistic numerical results. An inadequate mesh introduces severe discretization errors and solver instability, embodying the principle that a compromised grid inevitably produces compromised solutions.

Historically, practitioners evaluated grid adequacy by conducting grid dependency studies, progressively refining the mesh from several million to tens of millions of cells until designated physical quantities of interest asymptotically ceased to vary. While still encountered in academic literature, global grid refinement studies require substantial computational overhead and presume asymptotic convergence, an assumption that can fail in practice if solvers encounter oscillatory convergence behaviors.

A more targeted approach relies on physics-based evaluation criteria derived from the underlying flow models. In turbulent flow regimes, computing integral length scales highlights spatial regions dominated by interacting vortex structures where localized cell concentration is physically necessary. Furthermore, applying a grid refinement index based on local turbulent kinetic energy dissipation allows users to evaluate whether a mesh is sufficiently resolved. Standard engineering practices consider a grid physically well-resolved when this index value exceeds five, ensuring that elements are concentrated to capture physical gradients rather than arbitrarily reducing global truncation error.

// Core Geometric Mesh Quality Metrics

Regardless of the meshing software employed, several fundamental geometric metrics must be monitored to maintain solver stability and prevent unphysical numerical diffusion.

Non-orthogonality quantifies the angular deviation between the normal vector of a cell face and the vector connecting the centroids of the two adjoining cells sharing that face. As non-orthogonality increases, additional numerical diffusion is introduced into the diffusion terms of the discretized equations. If this deviation reaches ninety degrees, the mathematical volume associated with the face calculation collapses to zero, creating a singularity in the solver algorithms. Angles exceeding ninety degrees result in catastrophic negative cell volumes. In practical OpenFOAM workflows, a maximum non-orthogonality limit of eighty degrees is typically enforced. While a grid exhibiting eighty degrees of non-orthogonality contains visible geometric distortion, the numerical discretizations can reliably maintain convergence.

Face skewness evaluates the spatial distance between the geometric face center and the point where the vector connecting adjacent cell centers intersects that face. Skewed cells introduce artificial diffusion and non-physical oscillations into the solution fields. High skewness values frequently manifest near sharp geometric corners and aerofoil trailing edges. Although finite volume formulations incorporate explicit numerical skewness corrections, the standard quality limit in OpenFOAM restricts skewness to a maximum value of eight, ensuring that the face intersection point does not deviate excessively outside cell boundaries.

The cell aspect ratio represents the proportion between the longest and shortest characteristic dimensions of a three-dimensional cell. While large aspect ratios introduce directional numerical diffusion, highly anisotropic aspect ratios ranging from ten thousand to one hundred thousand are standard and necessary within boundary layer grids to capture steep normal gradients efficiently. However, to evaluate forces and skin friction with maximum accuracy, keeping the aspect ratio near the boundary as low as computationally feasible remains recommended.

Smoothness governs the spatial expansion rate between neighboring cells across refinement transitions. Abrupt size variations between adjacent elements generate significant local truncation errors. Standard best practice dictates that cell size changes between adjacent cells should not exceed twenty percent, which corresponds to a maximum expansion factor of 1.2.

// Cell Topologies and Meshing Characteristics

The geometric topology and spatial orientation of cells directly influence numerical accuracy. Aligning cell faces parallel to the primary flow direction minimizes numerical truncation errors, an alignment readily achieved using hexahedral or polyhedral elements. While historical computational limitations favoured structured hexahedra over tetrahedral elements, modern gradient schemes ensure that hexahedral, tetrahedral, and polyhedral meshes yield comparable results when numerical schemes are appropriately calibrated.

Within the native OpenFOAM framework, `snappyHexMesh` operates strictly as a hexahedral-dominant Cartesian meshing utility. It cannot natively construct unstructured tetrahedral or polyhedral meshes. Because octree-based refinement in `snappyHexMesh` imposes a strict two-to-one cell division ratio, the transition zones between differing refinement levels can introduce localized numerical dissipation. Managing these transitions using adequate buffer layers in `nCellsBetweenLevels` mitigates sudden geometric jumps across the computational domain.

// Automated Quality Enforcement in snappyHexMeshDict

Geometric mesh quality is continually monitored and enforced throughout the castellation, snapping, and layer addition stages by parameters defined within the `meshQualityControls` sub-dictionary of `snappyHexMeshDict`. Keywords such as `maxNonOrtho`, `maxInternalSkewness`, and `maxBoundarySkewness` establish strict geometric thresholds. Additional parameters, including `minVol`, `minDeterminant`, `minFlatness`, and `minTwist`, prevent the creation of collapsed, concave, or inverted cell geometries.

During the snapping and layer extrusion phases, any vertex displacement that causes local cells to violate these defined quality limits is systematically scaled back by the `errorReduction` factor over several smoothing iterations governed by `nSmoothScale`. If boundary layer inflation struggles to complete due to strict geometric thresholds, the nested `relaxed` sub-dictionary within `meshQualityControls` allows users to specify relaxed criteria that take effect after a designated number of iterations, balancing near-wall layer coverage with overall numerical validity.

// Mesh Validation with checkMesh and Repair Workflows

Following the completion of the meshing process, the structural and geometric integrity of the grid must be formally evaluated using the OpenFOAM `checkMesh` command-line utility. Executing `checkMesh` parses the complete mesh topology located in the `constant/polyMesh` directory and outputs a comprehensive statistical report detailing non-orthogonality, skewness, aspect ratios, and boundary patch counts.

If `checkMesh` reports a failed status for a specific metric, the solver is not mechanically prevented from attempting a run. However, a failure serves as a critical warning that numerical divergence or significant inaccuracies will occur unless discretization schemes are adjusted or the underlying mesh is improved. For cases where failure is caused by a few isolated, highly distorted cells at complex topological intersections, users can isolate their spatial coordinates and employ cell removal manipulation utilities to excise the defective elements. Nevertheless, manual cell removal represents an emergency intervention; achieving a high-quality simulation ultimately requires correcting underlying CAD surface defects, refining local background grids, and optimizing parameter settings in `snappyHexMeshDict`.

// EOF
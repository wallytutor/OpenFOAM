= Geometry <chapter-2>

In computational fluid dynamics (CFD), defining boundary conditions is not merely a software requirement but a core mathematical formulation. A CFD simulation is mathematically structured as an _initial boundary value problem (IBVP)_, meaning that to obtain a unique solution, the user must specify both initial conditions throughout the domain and precise boundary conditions at all physical boundaries. Because boundaries act as the primary physical drivers of the simulation—dictating exactly what enters, exits, or interacts with the domain—they must be highly realistic. Setting unphysical boundary values is one of the most common causes of solver instability and divergence. The first step towards a proper problem formulation is to define the geometry of the domain, the main subject of this chapter.

Before initiating the spatial discretization workflow in `snappyHexMesh`, a clean and mathematically well-defined representation of the target geometry is required. Although meshing guides do not focus primarily on fundamental computer-aided design principles, selecting the appropriate tool for geometry conception is a critical prerequisite. Modern engineering workflows can utilize various CAD suites depending on operational requirements, licensing constraints, and platform availability. Users seeking local and open-source software can employ tools such as FreeCAD or Salome to construct parametric geometries and perform topological operations. For those more prone to use code directly, geometric conception under Gmsh is also possible. Commercial CAD software suites are equally viable, as all standard engineering design packages serve the common primary purpose of defining the boundary representations of the physical model. Conversely, artistic modeling packages such as Blender are generally not recommended for engineering workflows, as they are not specifically tailored for the high precision and parametric constraints required by computational fluid dynamics configurations.

== Surface export formats

To be ingested by the native OpenFOAM meshing ecosystem, surface models must be exported from CAD software in triangulated formats, specifically Stereolithography (`.stl`) or Wavefront Object (`.obj`) representations. These files describe the outer boundaries of the physical domain using planar triangular facets. Regardless of the chosen export format, the geometric surface must be watertight, meaning it must form a completely closed manifold without unstitched edges, holes, or overlapping surfaces. In addition, all triangular surface face normals must be consistently and correctly oriented outward or inward relative to the active fluid domain to avoid ambiguity during cell classification. In the standardized case directory layout of modern OpenFOAM distributions, all surface geometry files must be placed within the `constant/geometry/` directory, which replaces the legacy `triSurface/` directory used in older versions of the software.

== Feature extraction and edge definition

Capturing sharp geometric features and distinct topological edges is essential for accurate surface snapping and boundary alignment. In OpenFOAM, this process is managed by an auxiliary pre-processing utility designated as `surfaceFeatures`. This utility analyzes the triangulated surface file according to parameters defined in the `surfaceFeaturesDict` configuration file. By evaluating user-specified inclusion angles, `surfaceFeatures` identifies curvature boundaries, intersections, and sharp angles across the surface geometry, writing the resulting edge data to an edge mesh file with an `.eMesh` extension inside the `constant/geometry/` directory. In addition to sharp edges, the utility can compute surface proximity and point closeness metrics, which provide geometric span data that can later be referenced within `snappyHexMeshDict` to guarantee adequate spatial resolution across narrow flow passages.

== Specifying geometric surfaces <ch2-geometry>

Within the master configuration dictionary `snappyHexMeshDict`, all surface geometries and spatial entities must be declared within the top-level `geometry` sub-dictionary. External CAD files are incorporated by declaring a user-defined name for the entity and assigning it the `triSurfaceMesh` type, followed by the specific file name of the imported STL or OBJ file. If the imported surface file contains distinct topological regions or grouped faces created during the CAD export phase, these can be explicitly mapped using a nested `regions` sub-dictionary. Mapping named regions allows the user to reassign individual geometric face groups to dedicated patch names, facilitating localized surface refinement and simplifying the subsequent application of specific numerical boundary conditions in the simulation setup. A minimal example is provided in @lst-geometry-cylinder.

#figure(
  ```cpp
  geometry
  {
      cylinder
      {
        type triSurfaceMesh;
        file "cylinder.stl";
      }
  }
  ```,
  caption: [Add file `constant/geometry/cylinder.stl` referring to patch name `cylinder` (used in meshing and boundary conditions).],
) <lst-geometry-cylinder>

== Analytical and searchable geometric entities

In addition to importing external triangulated CAD surfaces, the `geometry` sub-dictionary permits the definition of native analytical geometric shapes. These internal entities are defined directly using mathematical primitives, eliminating the requirement to export secondary CAD files for spatial control. Users can configure shapes #footnote[See #link("https://cpp.openfoam.org/v13/classFoam_1_1searchableSurface.html")[`searchableSurface`] for details.] such as `searchableBox`, `searchableSphere`, or `searchableCylinder` by specifying their canonical spatial coordinates, such as minimum and maximum diagonal vectors for bounding boxes, or center coordinates and radii for spherical volumes. The creation of some entities is illustrated in @lst-searchable-surfaces-example.

#figure(
  ```cpp
  geometry
  {
      // ... other surfaces here ...

      myBox
      {
          type      box;
          min       (0.0  0.0  0.0);
          max       (10.0 10.0 10.0);
      }

      mySphere
      {
          type       sphere;
          centre     (10.0  10.0  10.0);
          radius     5.0;
      }

      myCylinder
      {
          type      cylinder;
          point1    (0.0  0.0  0.0);
          point2    (0.0  0.0  1.0);
          radius    0.1;
      }
  }
  ```,
  caption: [Examples of searchable surfaces creation.],
) <lst-searchable-surfaces-example>

These searchable entities play multiple roles across the meshing workflow. Primarily, they are referenced within the downstream `castellatedMeshControls` section (discussed in @chapter-3) to establish volumetric refinement zones, enforcing higher cell densities in wakes, mixing layers, or designated flow paths. Furthermore, analytical geometries can serve as bounding templates or help define spatial domains when identifying simulation boundaries, ensuring that geometric structures and refinement zones are accurately linked throughout the execution of `snappyHexMesh`.

// EOF

// --------------------------------------------------------------------
// Laval nozzle geometry (as provided by Paul Bartkiewitz)
//
// Author: Walter Dal'Maz Silva
// Date: Mar 08 2025
// --------------------------------------------------------------------

// SetFactory("OpenCASCADE");

General.BackgroundGradient = 0;
General.Color.Background = {200, 220, 240};
General.Color.Foreground = White;
General.Color.Text = White;
General.Color.Axes = White;
General.Color.SmallAxes = White;
General.Axes = 0;
General.SmallAxes = 1;
Geometry.OldNewReg = 0;

// --------------------------------------------------------------------
// Reference coordinates / dimensions
// --------------------------------------------------------------------

Include "laval-points.geo";

// --------------------------------------------------------------------
// Lines
// --------------------------------------------------------------------

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};

Line(5) = {6, 7};
Line(6) = {7, 8};
Circle(7) = {8, 11, 9};
Line(8) = {9, 10};

Line(9) = {1, 6};
Line(10) = {5, 10};

// --------------------------------------------------------------------
// Surface
// --------------------------------------------------------------------

Line Loop(1) = {1, 2, 3, 4, 10, -8, -7, -6, -5, -9};
Plane Surface(1) = {1};

// --------------------------------------------------------------------
// Discretization
// --------------------------------------------------------------------

Transfinite Curve {1, -5} = 50;
Transfinite Curve {2, -6} = 20;
Transfinite Curve {3, -7} = 20;
Transfinite Curve {4, -8} = 100;
Transfinite Curve {-9, 10} = 50;

// --------------------------------------------------------------------
// Boundary layer
// --------------------------------------------------------------------

Field[1] = BoundaryLayer;
Field[1].PointsList = {6:10};
Field[1].CurvesList = {5:8};
Field[1].AnisoMax = 1.0e+10;
Field[1].Quads = 1;
Field[1].Ratio = 1.1;
Field[1].Size = 0.0001;
Field[1].SizeFar = 0.002;
Field[1].Thickness = 0.003;
BoundaryLayer Field = 1;

// --------------------------------------------------------------------
// Make symmetric wedge
// --------------------------------------------------------------------

// theta = 10 * Pi / 180;

// Rotate {{1, 0, 0}, {0, 0, 0}, theta} { Surface{1}; }

// ext[] = Extrude {{1, 0, 0}, {0, 0, 0}, -2*theta} {
//   Surface{1}; Layers{1}; Recombine; };

// Physical Volume("volume") = { ext[1] };
// Physical Surface("inlet") = {8};
// Physical Surface("outlet") = {3};
// Physical Surface("wall") = {4:7};
// Physical Surface("front") = {1};
// Physical Surface("back") = {9};

// --------------------------------------------------------------------
// Meshing
// --------------------------------------------------------------------

Mesh.Algorithm = 8;

Mesh.Smoothing = 1;
Mesh.RecombineAll = 1;
Mesh.MeshSizeMax = 0.001;

Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 1;
Mesh 2;

Save "laval-2d.unv";

// --------------------------------------------------------------------
// EOF
// --------------------------------------------------------------------

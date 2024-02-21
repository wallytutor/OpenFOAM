// -----------------------------------------------------------------------------
//
//  Geometry assembly
//
// -----------------------------------------------------------------------------

General.AbortOnError = 1;
General.BackgroundGradient = 0;
General.Color.Background = {0, 0, 0};
General.Color.Foreground = White;
General.Color.Text = White;
General.Color.Axes = White;
General.Color.SmallAxes = White;
General.Axes = 0;
General.SmallAxes = 1;
Geometry.OldNewReg = 0;
Geometry.Surfaces = 1;

CHARACTERISTIC = 0.039;

// -----------------------------------------------------------------------------
//
//  Shared parameters for project
//
// -----------------------------------------------------------------------------

// Box width [m].
L = 2.0;

// Box depth [m]
D = 4.0;

// Coordinate of outlet first point on the left [m].
d = 2.0;

//  Hypothetical material accumulated on bottom [m].
bedh = 0.050;

// Height of lower part of the box [m]. 
hl = 0.75;

// Height of upper part of the box [m]. 
ht = 2.05;

// Diameters of pipes [m].
d0 = 0.75;  // Air feed
d1 = 0.21;  // First bottom extraction
d2 = 0.27;  // Second bottom extraction

// Equivalent length of pipe openings in 2D [m].
e0 = Pi*d0^2/(4*L);
e1 = Max(Pi*d1^2/(4*L), CHARACTERISTIC);
e2 = Max(Pi*d2^2/(4*L), CHARACTERISTIC);

Printf("Opening 0 .... %.4f", e0);
Printf("Opening 1 .... %.4f", e1);
Printf("Opening 2 .... %.4f", e2);

// Plane perpendicular to air feed attack angle [rad].
theta = Atan(1198/720);
cost = Cos(theta);
sint = Sin(theta);

// Distance below feed after correcting diameter to equivalent length [m].
df = 0.08 + 0.5*(d0-e0);

// Thickness of boundary layer [m].
thickbl = 0.0800;

// First cell thickness over shell [m].
shellbl = 0.0025;

// -----------------------------------------------------------------------------
//
//  Geometry of outer shell of system
//
// -----------------------------------------------------------------------------

// Base line with extraction points, origin in the lower-left corner of box.
// Here we add an hypothetical `bedh` layer of accumulated material over the
// bottom wall so that traps that become too narrow in 2D idealization can
// actually capture the material. The angle of fall is 45° established by
// simply adding a point displaced of `bedh` wrt the trap holes.
p01 = newp; Point(p01)  = {0,                       bedh, 0};
p02 = newp; Point(p02)  = {1.038 - 0.5*e1 - 1*bedh, bedh, 0};
p03 = newp; Point(p03)  = {1.038 - 0.5*e1 + 0*bedh, 0,    0};
p04 = newp; Point(p04)  = {1.038 + 0.5*e1 + 0*bedh, 0,    0};
p05 = newp; Point(p05)  = {1.038 + 0.5*e1 + 1*bedh, bedh, 0};
p06 = newp; Point(p06)  = {3.010 - 0.5*e2 - 1*bedh, bedh, 0};
p07 = newp; Point(p07)  = {3.010 - 0.5*e2 + 0*bedh, 0,    0};
p08 = newp; Point(p08)  = {3.010 + 0.5*e2 + 0*bedh, 0,    0};
p09 = newp; Point(p09)  = {3.010 + 0.5*e2 + 1*bedh, bedh, 0};

// Last point of baseline is computed from angle and depth.
p10 = newp; Point(p10)  = {D - (hl-bedh)/Tan(42*Pi/180), bedh, 0};

// Wall at box depth `D`.
p11 = newp; Point(p11) = {D,      hl,          0};
p12 = newp; Point(p12) = {D,      hl+ht,       0};

// Outlet, outlet imbrication, and top wall.
p13 = newp; Point(p13) = {d,      hl+ht,       0};
p14 = newp; Point(p14) = {d,      hl+ht-0.052, 0};
p15 = newp; Point(p15) = {d-1.28, hl+1.198,    0};

// Particle feed requires some trigonometry... I hope not to need to recompute
// this otherwise a draft will be required to get the angles...
dy12 = df*sint + 0.5*cost + e0*sint;
p16 = newp; Point(p16) = {dy12/Tan(theta), hl+dy12, 0};
p17 = newp; Point(p17) = {(df+e0)*cost-0.5*sint, hl+(df+e0)*sint + 0.5*cost, 0};
p18 = newp; Point(p18) = {df*cost-0.5*sint, hl+df*sint+0.5*cost, 0};
p19 = newp; Point(p19) = {df*cost, hl+df*sint, 0};

// Closure point just above origin.
p20 = newp; Point(p20) = {0, 0.75, 0};

// Join all points.
l01 = newl; Line(l01) = {p01, p02};
l02 = newl; Line(l02) = {p02, p03};
l03 = newl; Line(l03) = {p03, p04};
l04 = newl; Line(l04) = {p04, p05};
l05 = newl; Line(l05) = {p05, p06};
l06 = newl; Line(l06) = {p06, p07};
l07 = newl; Line(l07) = {p07, p08};
l08 = newl; Line(l08) = {p08, p09};
l09 = newl; Line(l09) = {p09, p10};
l10 = newl; Line(l10) = {p10, p11};
l11 = newl; Line(l11) = {p11, p12};
l12 = newl; Line(l12) = {p12, p13};
l13 = newl; Line(l13) = {p13, p14};
l14 = newl; Line(l14) = {p14, p15};
l15 = newl; Line(l15) = {p15, p16};
l16 = newl; Line(l16) = {p16, p17};
l17 = newl; Line(l17) = {p17, p18};
l18 = newl; Line(l18) = {p18, p19};
l19 = newl; Line(l19) = {p19, p20};
l20 = newl; Line(l20) = {p20, p01};

// Distinguish shell from the rest.
Color Orange{ Curve{ l01, l02, l03, l04, l05, l06, l07, l08, l09, l10,
                     l11, l12, l13, l14, l15, l16, l17, l18, l19, l20 }; }

// Enclose region with a curve loop.
theloops[0] = newll;
Curve Loop(theloops[0]) = { l01, l02, l03, l04, l05, l06, l07, l08, l09, l10,
                            l11, l12, l13, l14, l15, l16, l17, l18, l19, l20 };

// -----------------------------------------------------------------------------
//
//  CALCULATION FOR INJECTION VELOCITY
//
// -----------------------------------------------------------------------------

a[] = Point{p19};
b[] = Point{p18};
dy = a[1] - b[1];
dx = a[0] - b[0];
slopeinlet = Atan(dy / dx);
Printf("Inlet angle = %.8f", slopeinlet);

// -----------------------------------------------------------------------------
//
//  INTERNAL FLAPS
//
// -----------------------------------------------------------------------------

Include "../make-flap.geo";

pars[] = {1.1, 1.0, 0.01, 1.2, 50.0};
Call AddInternalFlap;

pars[] = {3.1, 1.0, 0.01, 1.3, 15.0};
Call AddInternalFlap;

pars[] = {1.95, 1.9, 0.01, 1.2, 10.0};
Call AddInternalFlap;

// -----------------------------------------------------------------------------
//
//  Some discretization
//
// -----------------------------------------------------------------------------

ps1 = news; Plane Surface(ps1) = {theloops[]};

// Bottom traps.
Transfinite Curve{3, 7} = 5;

// Sides of bottom traps.
Transfinite Curve{2, -4, 6, -8} = 10;

// Inlet.
Transfinite Curve{17} = 15;

// Sides of inlet.
Transfinite Curve{16, -18} = 20;

// Outlet to cyclone. 
// Transfinite Curve{12} = 40;

// Some boundary layers (ordered from inlet).
// Transfinite Curve{19} = 10;
// Transfinite Curve{20} = 20;
// Transfinite Curve{1}  = 60;
// Transfinite Curve{5}  = 80;
// Transfinite Curve{9}  = 8;
// Transfinite Curve{10} = 50;
// Transfinite Curve{11} = 100;
// Transfinite Curve{13} = 5;
// Transfinite Curve{14} = 60;
// Transfinite Curve{15} = 40;

// Recombine.
// Recombine Surface {ps1};

// -----------------------------------------------------------------------------
//
//  Boundary modeling
//
// -----------------------------------------------------------------------------

// Top surface.
no = NEXTFIELD++;
Field[no] = BoundaryLayer;
Field[no].CurvesList = {13, 14, 15, 16};
Field[no].PointsList = {13, 17};
Field[no].Quads = 1;
Field[no].Ratio = 1.2;
Field[no].Size = shellbl;
Field[no].Thickness = thickbl;
BoundaryLayer Field = no;

// From entry to first drop.
no = NEXTFIELD++;
Field[no] = BoundaryLayer;
Field[no].CurvesList = {18, 19, 20, 1, 2};
Field[no].PointsList = {18, 3};
Field[no].Quads = 1;
Field[no].Ratio = 1.2;
Field[no].Size = shellbl;
Field[no].Thickness = thickbl;
BoundaryLayer Field = no;

// Bottom surface.
no = NEXTFIELD++;
Field[no] = BoundaryLayer;
Field[no].CurvesList = {4, 5, 6};
Field[no].PointsList = {4, 7};
Field[no].Quads = 1;
Field[no].Ratio = 1.2;
Field[no].Size = shellbl;
Field[no].Thickness = thickbl;
BoundaryLayer Field = no;

// Back surface.
no = NEXTFIELD++;
Field[no] = BoundaryLayer;
Field[no].CurvesList = {8, 9, 10, 11};
Field[no].PointsList = {8, 12};
Field[no].Quads = 1;
Field[no].Ratio = 1.2;
Field[no].Size = shellbl;
Field[no].Thickness = thickbl;
BoundaryLayer Field = no;

// -----------------------------------------------------------------------------
//
//  Preparing for OpenFOAM
//
// -----------------------------------------------------------------------------

Extrude {0, 0, 0.1} { Surface{ps1}; Layers{1}; Recombine; }

Physical Volume("volume") = {1};
Physical Surface("frontAndBack") = {1, 35};
Physical Surface("inlet") = {19};
Physical Surface("outlet-screw-1") = {5};
Physical Surface("outlet-screw-2") = {9};
Physical Surface("outlet-cyclone") = {14};
Physical Surface("walls") = {
    15, 16, 17, 18,      // Top
    20, 21, 22, 3, 4,    // Below inlet
    6, 7, 8,             // Bottom
    10, 11, 12, 13       // Back
};
Physical Surface("wing-1") = {23, 24, 25, 26};
Physical Surface("wing-2") = {27, 28, 29, 30};
Physical Surface("wing-3") = {31, 32, 33, 34};

// -----------------------------------------------------------------------------
//
//  Meshing controls
//
// -----------------------------------------------------------------------------

Mesh.SaveAll = 0;
Mesh.MshFileVersion = 2.2;

Mesh.Algorithm = 8;
Mesh.MeshSizeMax = CHARACTERISTIC;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 1;
Mesh.Smoothing = 20;

// -----------------------------------------------------------------------------
//
//  EOF
//
// -----------------------------------------------------------------------------

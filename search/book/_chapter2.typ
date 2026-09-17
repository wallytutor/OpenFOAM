#import "@preview/glossarium:0.5.10": make-glossary, gls
#show: make-glossary

= Fluid Dynamics <chapter-2>

#quote(
  block: true,
  attribution: [The Stone Roses, _Waterfall_ (1989).]
)["Stands on shifting sands, the scales held in her hands. The wind it just whips her and wails and fills up her brigantine sails."]

This chapter introduces fluid dynamics for CFD. It describes: governing equations, _i.e._ conservation of mass, momentum and energy; and, associated physical models, _e.g._ for viscosity, heat conduction and thermodynamics.

The equations describe fluid motion, forces and heat in time and three-dimensional (3D) space. Vector notation provides a mathematical framework to present the equations in a compact form. It enables the equations to be presented independently of any co-ordinate system, _e.g._ Cartesian (#gls("x")/#gls("y")/#gls("z")) or spherical (#gls("r")/#gls("theta")/#gls("phi")). It includes a standard set of algebraic operations, _e.g._ the inner (dot) and outer products.

The notation helps to ensure that the terms in equations are unchanged, or _invariant_, under a co-ordinate system transformation. Without invariance, a flow solution, _e.g._ along a pipe, would be dependent on the orientation of the pipe with respect to the co-ordinate system. Logically this dependence cannot exist; the laws of motion are the same in all "inertial frames".#footnote[Galileo Galilei, _Dialogo sopra i due massimi sistemi del mondo_, 1632.]

The derivation of the governing equations uses a control volume #gls("V") bounded by a surface #gls("S"), presented using the two-dimensional (2D) illustration above. We use $d #gls("V")$ and $d #gls("S")$ to describe an infinitesimally small volume and surface, respectively, and $T B D$ is the unit normal vector for each increment of surface $T B D$ , discussed in Sec. 2.1. It is important to note in any derivation whether the volume is defined as fixed in space or moving with the fluid.

Each derivation generally begins with a definite integral of some quantity, _e.g._ $T B D$ , over the volume $T B D$ denoted by

$ T B D $
If this notation is unfamiliar, understand it to mean a summation for all increments of volume $T B D$ that make up the total volume $T B D$ . The summed values are $T B D$ , where $T B D$ is the value in the respective $T B D$ .
The derivations also use integrals over the surface $T B D$ , _e.g._

$ T B D $
where $T B D$ . Volume and surface integrals are connected through Gauss’s Theorem, introduced in Sec. 2.4.

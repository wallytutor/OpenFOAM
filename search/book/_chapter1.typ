#import "@preview/glossarium:0.5.10": make-glossary, gls
#show: make-glossary

= Introduction <chapter-1>

#quote(
  block: true,
  attribution: [The Jam, _Start!_ (1980).]
)["If we get through for two minutes only, it will be a start."]

Fluid dynamics is concerned with the motion of fluids (liquids and gases) and the forces on them. "Computational" refers to computation of the flow and forces using numerical analysis. A literal definition of _computational fluid dynamics_ might therefore be _the prediction of fluid_ _motion and forces by computation using numerical analysis_.

Today, the term computational fluid dynamics, usually abbreviated to CFD, is used to describe a broader range of calculations for a wide variety of scientific and engineering applications. In particular, it is commonly used for applications involving heat, including the following:

- engines, _e.g._ internal combustion engines, turbines;
- heat recovery, _e.g._ heat exchangers;
- thermal management, _e.g._ cooling, exhaust systems;
- thermal comfort, with heating, ventilation, and air conditioning.


Thermodynamics is an important consideration in many of these applications. It relates internal energy to temperature, which affects the flow of heat. Further sources of heat include thermal radiation and chemical reactions, in particular combustion. Heat transfer may involve conduction in solid materials, coupled with the fluid flow, known as conjugate heat transfer.

If we include the examples mentioned above, a modern definition of computational fluid dynamics would be *the prediction of fluid motion and forces by computation using* *numerical analysis, generally extended to include heat, thermodynamics, chemistry and* *solids*.

Numerical analysis provides many methods and algorithms that are suitable for CFD. The methods include finite volume, finite element and finite difference, which calculate the distributions of properties, _e.g._ pressure, velocity and temperature, over regions of space which are usually fixed. Alternative methods attribute properties to particles represented by points in space, whose motions are calculated.

Particle methods are often used to approximate small scale flow features such as liquid sprays, e.g. for cooling, coating, cleaning, agriculture, food production, fire suppression, emission reduction and fuel injection. Solid particles can also be simulated in applications such as filtration, erosion, and fluidised beds.

As a presentation of general principles, this book describes numerical methodsto solve problems in fluid dynamics, up to and including heat and some basic thermodynamics, without extending into thermal radiation, chemistry and solids.

The book is solely dedicated to the _finite volume method_, the chosen method for many decades in the most popular, general-purpose CFD codes, including OpenFOAM. It also does not describe particle methods.

== Solution overview

Let us imagine calculating fluid flow along a pipe with CFD. To perform the calculation first requires a description of the problem by:

- the domain occupied by the fluid, _i.e._ the internal region of the pipe;
- equations that represent the fluid behaviour, in terms of properties such as pressure #gls("pressure") and velocity #gls("velocity");
- conditions at the boundary of the fluid domain and _initially_ within the domain for the fluid properties.

This description is represented in CFD by the following:

- a computational mesh for the fluid;
- "discrete" equations and algorithms to calculate #gls("pressure") and velocity #gls("velocity") ;
- boundary and initial conditions for #gls("pressure") and velocity #gls("velocity") .

@chapter-2 introduces the governing equations and basic models for fluid motion, forces and heat. Turbulence, which commonly occurs in many flows, is introduced in @chapter-6 and its standard modelling is described in @chapter-7.

The finite volume method is presented in @chapter-3 to express equations in discrete form using its geometric representation of a computational mesh. The algorithms used to solve the matrix equations and to couple sets of equations are given in @chapter-5.

@chapter-4 describes boundary conditions first from a numerical perspective, _i.e._ how they modify the matrix equations to influence the solution. It then covers a range of conditions that represent the behaviour at boundaries in many problems.

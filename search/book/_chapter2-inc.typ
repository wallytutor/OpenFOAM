This chapter introduces fluid dynamics for CFD. It describes: governing equations, _i.e._ conservation of mass, momentum and energy; and, associated physical models, _e.g._ for viscosity, heat conduction and thermodynamics.

The equations describe fluid motion, forces and heat in time and three-dimensional (3D) space. Vector notation provides a mathematical framework to present the equations in a compact form. It enables the equations to be presented independently of any co-ordinate system, _e.g._ Cartesian ( $T B D$ / $T B D$ / $T B D$ ) or spherical ( $T B D$ / $T B D$ / $T B D$ ). It includes a standard set of algebraic operations, _e.g._ the inner (dot) and outer products.

The notation helps to ensure that the terms in equations are unchanged, or _invariant_, under a co-ordinate system transformation. Without invariance, a flow solution, _e.g._ along a pipe, would be dependent on the orientation of the pipe with respect to the co-ordinate system. Logically this dependence cannot exist; the laws of motion are the same in all “inertial frames”.#footnote[Galileo Galilei, _Dialogo sopra i due massimi sistemi del mondo_, 1632.]

The derivation of the governing equations uses a control volume $T B D$ bounded by a surface $T B D$ , presented using the two-dimensional (2D) illustration above. We use $T B D$ and $T B D$ to describe an infinitesimally small volume and surface, respectively, and $T B D$ is the unit normal vector for each increment of surface $T B D$ , discussed in Sec. 2.1. It is important to note in any derivation whether the volume is defined as fixed in space or moving with the fluid.

Each derivation generally begins with a definite integral of some quantity, _e.g._ $T B D$ , over the volume $T B D$ denoted by

$ T B D $
If this notation is unfamiliar, understand it to mean a summation for all increments of volume $T B D$ that make up the total volume $T B D$ . The summed values are $T B D$ , where $T B D$ is the value in the respective $T B D$ .
The derivations also use integrals over the surface $T B D$ , _e.g._

$ T B D $
where $T B D$ . Volume and surface integrals are connected through Gauss’s Theorem, introduced in Sec. 2.4.


== Pressure

The equations of fluid dynamics in CFD treat the fluid as a continuous medium, or _continuum_. It is continuous in the sense that we consider the fluid as having no “empty space” by ignoring its molecular nature. We assume it has properties that vary from point to point and are continuous throughout the solution domain, and whose derivatives are also continuous.

Pressure is an important property of a fluid, denoted by $T B D$ . It describes the amount of force per unit surface area which _acts on_ a surface, in the direction _perpendicular_ to the surface.

Pressure $T B D$ is a _scalar_ that produces a force _vector_ with direction normal to the surface. When pressure is applied to one side of a segment of surface with area $T B D$ , the force $T B D$ points _away_ from that side by

$ T B D $
where $T B D$ is the vector of _unit length_, normal to the surface $T B D$ . The term _vector_ denotes a geometric entity with magnitude and direction; a surface can be represented by a _surface area_ _vector_ $T B D$ of magnitude $T B D$ and direction $T B D$ .
While pressure exerts a force _on a surface_, the fluid experiences a force which is compressive in nature, assuming $T B D$ is positive.

Pressure is measured in SI units of pascal $T B D$ . From Sec. 2.13 onwards, however, we generally use $T B D$ to represent the _kinematic_ pressure, in units $T B D$ , obtained by dividing the true pressure by a constant density $T B D$ .

=== Scalar fields

The majority of properties, _e.g._ pressure, temperature, energy, density, volume, _etc._, can be represented by a single number, or _scalar_. A _scalar field_ describes a scalar property, _e.g._ pressure, which varies from point to point across some spatial domain.

Point locations can be defined in any co-ordinate system of axes, _e.g._ Cartesian ( $T B D$ , $T B D$ , $T B D$ ), and in any orientation. A scalar field is _invariant_, meaning the scalar values are the same irrespective of the co-ordinate system used.

In this book, space and fields will be described in a co-ordinate system with right-handed rectangular Cartesian axes. The axes are constructed by defining an origin $T B D$ from which three lines are drawn at right angles to each other, termed the $T B D$ , $T B D$ , $T B D$ axes. A right-handed set of axes requires that looking down the $T B D$ axis with $T B D$ nearest, an arc from the $T B D$ axis to the $T B D$ axis is in a clockwise sense.



== Velocity

Like force, velocity $T B D$ is a vector with direction and magnitude, with SI units of $T B D$ . Using the vector $T B D$ to denote the position of a particle of fluid, its velocity is

$ T B D $

=== Vector fields

While $T B D$ can be used to denote a single velocity with magnitude and direction, it can also denote a _vector field_ of velocity which varies from point to point across a spatial domain. A vector is represented by 3 numbers, relating to the co-ordinate system being used, _e.g._ $T B D$ , $T B D$ , $T B D$ , in the Cartesian system.

While the magnitude and direction of a vector is fixed, it is not invariant since the 3 _values_ depend on the co-ordinate system used. We represent a vector without reference to the co-ordinate system by bold text, _e.g._ “ $T B D$ ” (compared to scalar “ $T B D$ ”).

=== Basic vector algebra

Addition and subtraction of 2 vectors is performed by operating on respective components. Subtraction of two vectors $T B D$ and $T B D$ is performed by

$ T B D $
Multiplication of any vector $T B D$ by a scalar $T B D$ is performed by multiplying all the components by the scalar, _e.g._
$ T B D $
Addition and multiplication are _commutative_, _i.e._ variables can be in any order ( $T B D$ ). Subtraction is not commutative. Products between scalars and vectors are _distributive_, _i.e._
$ T B D $
and with additional scalar $T B D$ ,
$ T B D $
Division between a vector $T B D$ and a scalar is only relevant when the scalar is the second argument of the operation, _i.e._
$ T B D $



== Flow through a surface

The concept of flow through a surface appears in many areas of CFD, including fluid dynamics equations, numerical methods, boundary conditions and general flow calculations. When we talk about something that travels through a surface, the term _flux_ is generally used.#footnote[Joseph Fourier used the term _fluxion_ in relation to flow of heat in _Th__éorie analytique de la chaleur_, 1822; James Clerk Maxwell used the term _flux_ in _A Treatise on Electricity and Magnetism_, 1873.]

To quantify the flux of some property, we multiply the area of surface by the property at the surface. If the property is a vector, we take the component _normal to the surface_. For example, the flux $T B D$ associated with velocity through a surface segment of area $T B D$ would be $T B D$ .

As shown in the figure, $T B D$ can be calculated from $T B D$ by the inner product with normal vector $T B D$ of unit length, expressed as

$ T B D $
The flux associated with $T B D$ is
$ T B D $
It is a good habit to write $T B D$ first since the order matters with a tensor, _e.g._ stress $T B D$ , introduced in Sec. 2.6, _i.e._ $T B D$
=== Inner product of two vectors

The normal component of velocity $T B D$ is described in Eq. (2.5) by the inner, or “dot”, product of $T B D$ and $T B D$ . It is calculated for vectors $T B D$ and $T B D$ as shown in the figure below, where $T B D$ denotes the magnitude of the vector and $T B D$ is the internal angle between the two vectors.

The inner product of two vectors is a scalar invariant, since the magnitudes and angle are the same irrespective of the co-ordinate system used. It is calculated from vector components as follows:

$ T B D $

The inner product of two vectors is commutative, _i.e._ $T B D$ . It is distributive, _i.e._ with an additional vector $T B D$ ,

$ T B D $

Scalar multiplication and inner products are associative, _i.e._

$ T B D $

An inner product of a vector with itself is simply the square of the vector magnitude, _i.e._

$ T B D $



== Conservation of mass

The law of conservation of mass can be written

$ T B D $
where $T B D$ is the mass density of the fluid. The equation can be derived by considering a volume $T B D$ _fixed in space_ (note !!), bounded by a surface $T B D$ . If the volume is filled by a fluid with density $T B D$ , its mass is $T B D$
The rate of increase of mass inside the volume must match the rate of inflow of mass across the volume’s surface. The latter is calculated by integrating the mass flux $T B D$ over the surface, noting the negative sign due to $T B D$ pointing out of the volume.

Gauss’s Theorem relates surface and volume integrals by

$ T B D $
Equating the rate of mass increase to rate of inflow and applying Eq. (2.9) gives
$ T B D $
Since the integral is valid for any volume $T B D$ , it follows that the integrand (in $T B D$ ) must equal 0, resulting in Eq. (2.8).
=== Divergence

Divergence, denoted by $T B D$ , indicates the tendency of a vector field to point outward of a closed surface. For example, When the divergence of velocity $T B D$ is positive, the fluid is expanding; negative divergence indicates contraction. Imagine a volume of fluid with a complex distribution of $T B D$ at its bounding surface, below left. The $T B D$ calculation will isolate the diverging component (right) from uniform flow (centre).

The divergence of velocity is calculated by integrating — _i.e._ summing — fluxes over the closed surface. To define divergence _at a point_, we consider the limiting case where the volume tends to zero. For a surface $T B D$ that encloses a volume $T B D$ , divergence is the flux across the surface per unit volume, as $T B D$ , _e.g._

$ T B D $

=== The nabla operator

The nabla symbol $T B D$ can be considered a _vector operator_

$ T B D $
Divergence is the inner product with $T B D$ , _e.g._
$ T B D $



== Time derivatives

The conservation of mass Eq. (2.8) included a partial derivative in time $T B D$ relating to a _fixed region of space_. This is the _local rate of change of_ $T B D$ , relating to the change in $T B D$ in the fluid measured by an observer at a fixed location. It is not the time rate of change experienced by a mass of fluid particles as they move through space. In the same way, $T B D$ is not the acceleration experienced by the fluid.

Acceleration relates to the _material_, or _substantive_, derivative which describes the time rate of change of a _fixed mass of moving material_. It is denoted by $T B D$ and is related to the local rate of change, using $T B D$ as an example tensor of any rank, by

$ T B D $

The relation is derived from the _chain rule of differential_. In one dimension, it is illustrated by two particles of fluid that occupy positions 1 and 2 at some initial time, then positions 3 and 4 at a later time $T B D$ . The particles move in the $T B D$ -direction at speed $T B D$ , such that the particle at 1 later occupies the position 3 and the particle at 2 occupies position 4.

The material time derivative of $T B D$ , following the mass from 2 $T B D$ 4, is the sum of: the local change in $T B D$ at a fixed position $T B D$ (2 $T B D$ 3); and, the change due to the _gradient_ of  $T B D$ between positions 3 and 4, fixing time $T B D$ . This equates to

$ T B D $
The material derivative relation in Eq. (2.14) is simply the 3D equivalent of Eq. (2.15).
=== Gradient

The last term in Eq. (2.14) introduces the _gradient_ denoted by $T B D$ . If $T B D$ is a scalar, the gradient produces a vector whose magnitude and direction is that of the _steepest_ _gradient_.

The figure above illustrates the gradient using a surface that represents a distribution of a scalar field in 2 directions. The gradients at 3 locations are in the direction of steepest ascent.

When $T B D$ is a vector, the gradient produces a tensor, representing the direction and magnitude of steepest ascent for each of the 3 components of the vector.



== Forces at a surface

The next law to define is that of conservation of momentum, _i.e._ Newton’s second law of motion ( $T B D$ ) for fluids. It involves forces within the fluid so requires a description of forces at a surface $T B D$ bounding a volume $T B D$ .

The force $T B D$ is in the direction of the traction vector with magnitude of (traction $T B D$ surface area) — compare with Eq. (2.3)

$ T B D $
The tractions at a surface bounding some volume of fluid (or any continuum, _e.g._ solid) _depend_ _on the orientation of the surface_. Consequently, we cannot define forces within a fluid as traction vectors _at points_ within the fluid.
Instead, we require 3 _traction vectors_ $T B D$ , $T B D$ and $T B D$ , defined in planes perpendicular to one another, _i.e._ $T B D$ , $T B D$ and $T B D$ . This results in the _stress tensor_ $T B D$ with 9 components, consisting of 3 traction vectors, each containing 3 components.

The traction can be calculated at a surface with any orientation by taking the inner product of the unit normal vector $T B D$ and the stress tensor $T B D$ such that

$ T B D $

=== Tensors

We have defined the stress tensor#footnote[More precisely the “Cauchy stress tensor”, introduced by Augustin-Louis Cauchy in _De la pression ou_ _tension dans un corps solide_, 1827.], an entity with 9 component values, corresponding to our $T B D$ , $T B D$ and $T B D$ axes (or more specifically _base vectors of unit length_ aligned with our $T B D$ , $T B D$ and $T B D$ axes).

In fact, the term “tensor” describes any entity with multiple component values corresponding to the dimensions of space — here 3. A tensor has _rank_ $T B D$ , such that the number of component values for 3D space = $T B D$ .

In this book we use the term “tensor” to mean “tensor of rank 2” unless otherwise noted. A vector is a tensor of rank 1 and a scalar is rank 0.

The inner product of a vector $T B D$ and tensor $T B D$ produces a vector whose 3 components are evaluated as follows:

$ T B D $
This inner product tensor is commutative only if $T B D$ is symmetric since $T B D$ — see Eq. (2.24) for the _transpose_ $T B D$ .


== Conservation of momentum

The law of conservation of momentum can be written#footnote[Cauchy’s first law of motion (1827)]

$ T B D $
where $T B D$ represents any body force per unit mass. Body forces represent any force which does not act at a bounding surface, including those that act at a distance, such as gravitational force.
The equation is derived by considering the time rate of change of the momentum of a mass of particles. We consider a _volume of material of fixed mass_ moving through space and therefore present rate of change by the material derivative $T B D$ .

Applying Gauss’s theorem and Eq. (2.16), the surface force is

$ T B D $
Equating the rate of change of momentum to the forces and, noting mass is fixed so $T B D$ is constant in time, gives
$ T B D $
The integrand must equal 0, resulting in Eq. (2.19).


Divergence was described in section 2.4 as the flux across a surface per unit volume as $T B D$ . The divergence of stress similarly represents the stress flux across the surface, _i.e._ force, per unit volume as $T B D$ , given by

$ T B D $
Eq. (2.19) specifically relates to _linear_ momentum. Instead, conservation of _angular_ momentum#footnote[Cauchy’s second law of motion (1827)], in the absence of any “couple stresses” that generate a moment field, is given by
$ T B D $
The derivation of Eq. (2.23) is fairly complex so is omitted here.
=== Tensor symmetry

Symmetry in a tensor $T B D$ refers to components being symmetric about the diagonal, _i.e._ $T B D$ , $T B D$ and $T B D$ . The _transpose_ of a tensor, denoted by the ‘ $T B D$ ’ superscript, switches components across the diagonal such that:

$ T B D $
$T B D$ is therefore _symmetric_ if $T B D$ . A _skew_ (anti-symmetric) tensor has $T B D$ . A tensor can be decomposed into symmetric and skew parts by:
$ T B D $



== Flow in a volume

We presented conservation of momentum, Eq. (2.19), in terms of the material derivative $T B D$ . To solve the equation in this form, would require a method that tracks particles of fluid as they move around.

Generally it is easier to solve the equation by _fixing_ the volume of space and solving for the fluid motion through it. To enable this, we replace the material derivative for $T B D$ using the following expression:

$ T B D $
Eq. (2.26) is derived from the material derivative Eq. (2.14) by
Conservation of momentum can then be written

$ T B D $
A term of the form $T B D$ or $T B D$ represents the bulk motion of the fluid at velocity $T B D$ , which transports the property $T B D$ (a tensor of any rank) by advection#footnote[The term “convection” is sometimes used, but it has conflicting meanings, so we avoid it here.]. For example, if the property $T B D$ is temperature, then advection will transport heat from one region of the flow domain to another.
Advection contains a divergence derivative so represents a flux across a surface per unit volume. The example in Eq. (2.27) represents the flux of momentum, where the advected property is itself $T B D$ (or $T B D$ , depending whether $T B D$ is associated with bulk flow or the advected property).

=== Outer product of two vectors

The advection term in Eq. (2.27) includes a product of two vectors $T B D$ . This is the outer product of two vectors#footnote[sometimes written “ $T B D$ ”, but we write it as a product with no symbol, similar to a scalar multiplication.] which produces a tensor by

$ T B D $

=== Inner product of two tensors

The inner product of two tensors, _e.g._ $T B D$ produces a tensor $T B D$ where the components are (replacing $T B D$ with $T B D$ )

$ T B D $

=== Identity tensor

The identity tensor $T B D$ is the tensor equivalent of unity (one) such that for any tensor $T B D$ ,

$ T B D $



== Conservation and boundedness

The material derivative Eq. (2.14) includes the $T B D$ term. However, the term can take for the form $T B D$ , by applying conservation of mass Eq. (2.8), as shown in Sec. 2.8. The two terms are related by the product rule

$ T B D $

While the two terms $T B D$ and $T B D$ appear similar, they affect the solution of an equation is different ways.

The $T B D$ term _maintains conservation_ in $T B D$ . All variables are inside (to the right of) the divergence $T B D$ , so when we integrate it over a volume $T B D$ , it can be entirely transformed by Gauss’s theorem to an integral of the flux of $T B D$ at $T B D$ .

If we split $T B D$ into two sub-volumes $T B D$ and $T B D$ , the surface integrals are also split over two surfaces $T B D$ and $T B D$ , where $T B D$ is the part of the surface common to both volumes.

The fluxes at $T B D$ for each sub-volume have equal magnitude but opposite sign, since their respective surface area vectors $T B D$ point outwards from the volume in opposing directions. These fluxes cancel one another out, such that the sum of integrals over $T B D$ and $T B D$ is identically equal to the integral over $T B D$ .

Thus, the term $T B D$ ensures conservation in $T B D$ across a surface separating regions of the flow domain. Conservation is guaranteed at all points in the limit $T B D$ for any sub-volume.

The $T B D$ form cannot transform to a surface integral so does not ensure conservation. Instead it _ensures boundedness_, as demonstrated by the solution of $T B D$ in the following equation in one ( $T B D$ ) dimension (1D)

$ T B D $
The solution to this equation has the form $T B D$ . If $T B D$ is uniform, then the profile of $T B D$ does not change but simply translates at a speed $T B D$ as shown in the left diagram. If $T B D$ varies with $T B D$ , the profile changes, _e.g._ flattens, as shown in the right diagram.
In both cases, the $T B D$ profile _remains within the original bounds_ of $T B D$ and $T B D$ ; the solution is said to be _bounded_. While the behaviour is illustrated in 1D, it extends to 3D.

Thus, the $T B D$ term ensures boundedness. By contrast, the $T B D$ term, while ensuring conservation, does not ensure boundedness when $T B D$ . The two terms are connected by $T B D$ which changes $T B D$ due to expansion/contraction of the fluid, as discussed in Sec. 2.4. Its effect is illustrated by flattening of the profile in the right diagram, since a non-uniform $T B D$ corresponds to $T B D$ .



== Fluid deformation

Conservation of momentum Eq. (2.19) includes the $T B D$ term describing the forces within a fluid. At a macroscopic level, the forces exist due to interactions between fluid particles and change due to deformation of the fluid as it moves.

From the velocity field $T B D$ , we need to isolate pure deformation of the fluid from other characteristics of the flow. The figure below shows a rectangular fluid element under shear in the $T B D$ -direction due to a uniform velocity gradient in the $T B D$ -direction.

The velocity gradient $T B D$ calculates _some_ rate of shear of the rectangular element. The upper diagram shows shear in one direction, known as _parallel shear_.

Parallel shear can be decomposed into pure shear (left) and a rotation (right). The rotational component is represented by the skew part of $T B D$ , see Eq. (2.25). Since a rotation does not involve deformation, it must be removed from any measure of it, which leaves the symmetric part of the $T B D$ tensor.

Therefore, the rate of deformation tensor $T B D$ and spin (rotational) tensor $T B D$ are defined as

$ T B D $

=== Trace of a tensor

The trace of a tensor $T B D$ is the sum of the diagonal components, producing a scalar

$ T B D $
Where the tensor is $T B D$ , we conclude from Eq. (2.13) that
$ T B D $

=== Deviatoric and spherical tensors

A tensor can be decomposed into the sum of its deviatoric part and a spherical part as follows:

$ T B D $
where $T B D$ is the identity tensor, see Eq. (2.30). The deviatoric part subtracts the mean of the trace from each diagonal component such that the resulting tensor is “trace-free”, _i.e._ $T B D$ .


== Vorticity

Vorticity describes the tendency for a fluid to rotate locally, defined as

$ T B D $
where the $T B D$ operator is the _curl derivative_. The curl of a vector $T B D$ is evaluated in Cartesian co-ordinates by
$ T B D $
The local rotation which vorticity measures is not due to the flow following a curved path. In fact, vorticity is often associated with shear when streamlines may be straight and parallel, _i.e._ not curved at all.
Vorticity is difficult to picture under shear because the deformation masks the local rotation. By separating the deformation, as in the figure in Sec. 2.10, the local rotation is revealed.

Vorticity $T B D$ is often demonstrated by a vortex ring produced by an air “cannon” with smoke to visualise the flow. It reveals flow _circulation_ around sections of the torus, from front to back. The vorticity vectors are normal to the planes of circulation, along the axis of the torus.

Stokes’s theorem#footnote[posed by George Stokes as a examination question for the 1854 Smith’s Prize for progress in mathematics and natural philosophy at the University of Cambridge, awarded to James Clerk Maxwell and Edward Routh.] relates _circulation_ — the integral of $T B D$ around a closed curved line, $T B D$ — to the integral of vorticity $T B D$ over a section of surface $T B D$ bounded by the curve, according to:

$ T B D $
In Eq. (2.39), $T B D$ is a vector representing a segment of the line $T B D$ . As you stand on $T B D$ looking in the direction $T B D$ , with your head in the direction $T B D$ , $T B D$ is oriented to your left.
Vorticity is related to the spin tensor $T B D$ in Eq. (2.33) by $T B D$ , where $T B D$ is the Hodge dual operator which extracts components of a vector from a tensor $T B D$ as shown below:

$ T B D $


The term “spin tensor” emphasises _local_ rotation as opposed to the general path in a similar way that the spin on a ball is clearly distinguished from its trajectory, which curves due to gravity.


== Newtonian fluid

In Sec. 2.6 we introduced forces in a fluid. At the macroscopic scale of a continuum, we characterise a fluid’s response to applied forces through _constitutive models_.

The Newtonian (or linear viscous) fluid is the most common constitutive model that represents the behaviour of many liquids and gases. It states that a fluid at rest (or uniform velocity) does not sustain shear stress  $T B D$ ; it can be expressed by

$ T B D $

The model is a continuum representation#footnote[attributed to Adhémar Jean Claude Barré de Saint-Venant, 1843, and George Gabriel Stokes, 1845; derived previously using molecular models by Claude-Louis Navier, _M__émoire sur les lois du mouvement des_ _fluides_, 1822.] of Newton’s law of viscosity#footnote[Isaac Newton, _Philosophiae Naturalis Principia Mathematica_, 1687.] which states that shear stress $T B D$ is proportional to velocity gradient by the dynamic viscosity $T B D$ .

It is a specific case of the more general Stokesian fluid, defined as $T B D$ , where $T B D$ is deformation rate, Eq. (2.33). The Newtonian model assumes: (1) the fluid is isotropic, _i.e._ the value of $T B D$ is independent of the direction in which it is measured; (2) zero _bulk_ viscosity associated with a change in volume o the fluid.

Using $T B D$ in Eq. (2.41) ensures that the shear stress is induced by deformation only. Taking the deviatoric part, $T B D$ , ensures viscous stresses are not generated by volume changes, which are represented by $T B D$ .

This is due to the “ $T B D$ ” operator subtracting $T B D$ from _each_ diagonal component of $T B D$ , giving a total of $T B D$ ( $T B D$ ) from all 3 diagonal components. From Eq. (2.35), $T B D$

=== Pressure gradient

When we substitute Eq. (2.41) in $T B D$ in Eq. (2.19), the pressure part is $T B D$ . The term $T B D$ is equivalent to gradient $T B D$ .

Like divergence of stress in Eq. (2.22), the gradient of pressure represents the pressure flux across the surface per unit volume as $T B D$ , according to

$ T B D $
Since $T B D$ holds for any variable $T B D$ , a gradient term is conservative, like divergence, and can be converted to a surface integral under an equivalent Gauss’s Theorem
$ T B D $



== Incompressible flow

We can derive a set of equations for incompressible flow that includes:

- mass conservation, Eq. (2.8);
- momentum conservation, Eq. (2.19);
- the material derivative for $T B D$ , Eq. (2.26);
- the rate of deformation tensor, Eq. (2.33);
- the Newtonian fluid model, Eq. (2.41).


Combining Eq. (2.41) and Eq. (2.33) and using Eq. (2.36) for the deviatoric part of a tensor, gives an expression for stress:

$ T B D $
Substituting $T B D$ into Eq. (2.19) and applying Eq. (2.26) gives an equation for momentum for a Newtonian fluid:
$ T B D $
The derivation of the equation above uses Eq. (2.35) and the identity $T B D$ on page 84.
An incompressible fluid exhibits $T B D$  constant over time, _i.e._ for _moving volumes of fluid_ $T B D$ Combining the material derivative Eq. (2.14) and mass conservation Eq. (2.8) gives $T B D$ which results in the _incompressibility condition_

$ T B D $

A _homogeneous_, incompressible material exhibits $T B D$  = constant uniformly throughout the _entire fluid_. With that assumption, Eq. (2.45) can be written as

$ T B D $
This is the momentum equation for a homogeneous, incompressible Newtonian fluid. It includes: _kinematic_ viscosity $T B D$ in SI units of $T B D$ ; and $T B D$ represents _kinematic_ pressure, _i.e._ divided by $T B D$ , in SI units of $T B D$ .
The identity $T B D$ and the incompressibility Eq. (2.46) yield the terms in Eq. (2.47).

=== Pressure equation

Mass and momentum conservation, represented by Eq. (2.46) and Eq. (2.47) respectively, provide two equations — one scalar, one vector — for two fields, $T B D$ and $T B D$ . However, Eq. (2.46) cannot be solved in its own right since it provides only one equation for vector $T B D$ , containing 3 components.

A scalar equation, including both $T B D$ and $T B D$ , can be derived by taking _the divergence_ of Eq. (2.47) and eliminating terms by substituting Eq. (2.46), noting that $T B D$ . For $T B D$ and $T B D$ that are constant and uniform, the equation is

$ T B D $
This pressure equation can replace Eq. (2.46) to provide a pair of equations, with Eq. (2.47), for both $T B D$ and $T B D$ . Chapter 5 describes the algorithms used in the finite volume method which couple the two equations using a modified form of Eq. (2.48).


== Diffusion

The momentum equation for a homogeneous, incompressible, Newtonian fluid is presented in Eq. (2.47). In the case of zero body force, $T B D$ , and $T B D$  constant, Eq. (2.47) becomes

$ T B D $
This is a typical _transport equation_ characterised by:
- the local rate of change $T B D$ , described in Sec. 2.5;
- advection of $T B D$ by $T B D$ , described in Sec. 2.8;
- diffusion of $T B D$ by $T B D$ ;
- a “source” due to $T B D$ , described in Sec. 2.12.


Note that, since $T B D$  constant, $T B D$ , where $T B D$ denotes the Laplace operator.

The $T B D$ term is a special form of divergence in which the flux includes the _surface_ _normal gradient_ denoted by the operator $T B D$ where $T B D$ .

The term models _diffusion_ across the surface of the fluid element. Diffusion generally represents the transport of a fluid property — here, momentum — due to _fluctuating motions_ that are not captured by the bulk motion that is represented by the _continuum_ velocity $T B D$ .

Fluctuations include any random motion of particle constituents of matter, _e.g._ molecules, and turbulent structures. Through these motions, particles can pass across a surface boundary, transporting property $T B D$ through a gradient of $T B D$ (above, left).

Particles carrying higher $T B D$ move into regions of particles with lower $T B D$ and vice versa. Through particle collisions, high values of $T B D$ tend to reduce and low values increase (right).

=== Laplacian

“Laplacian”#footnote[After Pierre-Simon Laplace, _Th__éorie des attractions des sph__éro__ïdes et de la figure des_ _plan__ètes_, 1785 but described in 1761 both by Leonhard Euler and Jean-Baptiste le Rond d’Alembert.] describes a term of the form $T B D$ where $T B D$ is a diffusivity coefficient.#footnote[More precisely, “Laplacian” denotes a term without $T B D$ , _i.e._ $T B D$ , represented as $T B D$ but we include $T B D$ here because it is generally needed.] A Laplacian term is _conservative_ since all variables are to the right of a divergence, as described in Sec. 2.9. It is also _bounded_ since it tends to decrease high values and increase low values as shown above.

The Laplacian represents a flux due to $T B D$ across the surface, per unit volume, as $T B D$

$ T B D $



== Conservation of energy

The law of conservation of energy#footnote[from Hermann von Helmholtz, _Über die Erhaltung der Kraft_, 1847, who credited it equally to Julius von Mayer, 1842, and James Joule, 1843. “Joule’s apparatus” is the best know demonstration in which the work of a falling weight is converted into heat in water. The first law of thermodynamics is the law of conservation of energy applied to a system.] states that the total energy of an isolated system remains constant over time. Energy is not created or destroyed, but is transformed from one form to another. If we consider mechanical and thermal energy, it can be expressed by

$ T B D $
The equation describes the rate of change of mechanical and thermal energies for a volume of material of _fixed mass_ of fluid particles. Thermal energy is described by the specific internal energy $T B D$ , where _specific_ denotes _per unit mass_. Mechanical energy is represented by specific kinetic energy $T B D$ .
The rate of change of these combined energies is due to the input of mechanical power and heat from both surface fluxes and internal sources.

Mechanical power is calculated by the inner product of force and $T B D$ . The power over a surface segment $T B D$ is $T B D$ , representing a change in strain energy. The internal power source is $T B D$ , representing a change in potential energy due to general body forces.

Heat is provided by an internal source of strength $T B D$ per unit mass and a surface heat $T B D$ per unit surface area. The surface heat flux relating to heat _input_ takes a negative sign, $T B D$ , since $T B D$ follows the convention of pointing out of the volume.

Equation 2.51 is derived by equating integrals of rate of change of energy to the mechanical power and heat input and taking $T B D$ inside the volume integral since $T B D$ is constant in time, as follows:

$ T B D $
Gauss’s theorem converts the surface integrals of Eq. (2.52) into the following volume integrals:
$ T B D $
Substituting Eq. (2.53) for the surface integrals into Eq. (2.52) and equating the integrand to 0, as described previously in Sec. 2.7, yields Eq. (2.51).


== Temperature

In the conservation of energy Eq. (2.51), the mechanical kinetic energy, power flux and sources can be calculated from $T B D$ , $T B D$ and $T B D$ from the momentum Eq. (2.19). Heat sources can contribute to $T B D$ , _e.g._ from thermal radiation, chemical reactions _etc._

That leaves the heat flux term $T B D$ which represents _conduction_ of heat. It is commonly modelled by Fourier’s law#footnote[Joseph Fourier, _Th__éorie analytique de la chaleur_, 1822.] which states $T B D$ is proportional to the negative gradient of temperature $T B D$ , _i.e._

$ T B D $
where the constant of proportionality is thermal conductivity $T B D$ .
=== Temperature scale

The heat flux Eq. (2.54) requires temperature to be defined and measurable. Measurement requires a _scale_. Empirical scales correlate temperature with a measured physical property of a _working substance_, _e.g._ EMF at a junction of two metal alloys. Empirical scales have the drawbacks of: being dependent on the working substance; and, not actually defining temperature.

#footnote[in accordance with the zeroth law of thermodynamics.]

Instead, the _thermodynamic scale_ defines temperature as a measure of the average kinetic energy of random motions of particle constituents of matter. It provides an _absolute_ measure of temperature that is independent of the choice of working substance and includes a zero point#footnote[in accordance with the third law of thermodynamics.]. It must be measured in units with a zero point, such as the SI unit Kelvin, $T B D$ .

Substitution of our model Eq. (2.54) into Eq. (2.51) yields the term $T B D$ . It is logical that this is a Laplacian term since it represents diffusion which is associated with random motions of submicroscopic particles, as we we established in Sec. 2.14.

=== Ideal gas

The behaviour of many gases under typical working conditions is captured by the ideal gas equation of state

$ T B D $
where the specific gas constant is $T B D$ . It is calculated from the Universal Gas Constant $T B D$ in SI units, and the molar mass $T B D$ of the gas, with units of $T B D$ .
The ideal gas equation originates from classical thermodynamics as a combination of empirical laws#footnote[Benoît Clapeyron, _M__émoire sur la puissance motrice de la chaleur_, 1834.]. Later, it was derived from first principles from both statistical thermodynamics and kinetic theory, with temperature representing average kinetic energy.

The derivations assume that molecules have no volume, undergo purely elastic collisions and there are no inter-molecular forces.

A scale of temperature defined by the ideal gas equation of state is exactly equivalent to the thermodynamic temperature scale.



== Internal energy

The conservation of energy Eq. (2.51) introduces the internal energy using the specific quantity $T B D$ , measured in $T B D$ in SI units. It represents the total molecular energy consisting of: kinetic energy associated with temperature; potential energy due to particle forces, both _within_ particles as chemical bonds, and _between_ particles, _e.g._ van der Waals forces.

The energies at different scales can be summarised as:

- bulk — kinetic $T B D$ due to bulk motion, potential due to forces $T B D$ and $T B D$ ;
- molecular — kinetic characterised by $T B D$ , potential due to bonds.


To understand how energy is transferred between these scales, we can derive an equation for internal energy from Eq. (2.51), by cancelling terms in mechanical energy formed by taking the inner product of Eq. (2.19) with $T B D$ .

$ T B D $

The identity $T B D$ is the key element of the analysis. The $T B D$ term in Eq. (2.56) must represent the contribution of mechanical power to the internal energy, _i.e._ passing from bulk to molecular scale. Conversely, the $T B D$ term must contribute to mechanical energy.

Substituting $T B D$ from Eq. (2.41), Eq. (2.56) becomes

$ T B D $
The sign of $T B D$ depends on $T B D$ , _i.e._ whether the fluid is expanding or contracting. Since the sign can change, it therefore represents a _recoverable_ contribution to internal energy.
If we substitute the Newtonian model from Eq. (2.41)

$ T B D $
#footnote[This identity relies on the substitution of $T B D$ by $T B D$ which is equivalent due to symmetry of $T B D$ . $T B D$ can then be replaced by $T B D$ since $T B D$ .]The $T B D$ term is always positive since all components of $T B D$ are squared. Its contribution to $T B D$ is therefore _non-recoverable_ and represents mechanical power that is dissipated as heat. In the majority of CFD analyses, the $T B D$ contribution is small and can be ignored.
=== Double inner product of two tensors

The double inner product of two tensors, denoted by “ $T B D$ ”, is introduced in Eq. (2.56). It produces a scalar which is evaluated as the sum of the 9 products of the tensor components, eg:

$ T B D $
For a scalar $T B D$ , $T B D$ .


== Heat capacity

When the heat conduction model of Eq. (2.54) is substituted into Eq. (2.51), we produce an equation with two variables $T B D$ and $T B D$ . For example, starting with Eq. (2.56) for internal energy, and ignoring mechanical power and heat source terms, gives

$ T B D $
To be able to solve this equation, a relation between $T B D$ and $T B D$ is needed. _Heat capacity_ provides this, describing an amount of heat required to produce a unit change in $T B D$ . We use _specific_ heat capacity, relating to a constant volume process, whose SI units are $T B D$ , defined as
$ T B D $
Statistical mechanics provides a theorem known as _equipartition_ which provides quantitative predictions of $T B D$ .#footnote[originating from heat capacity studies by Alexis Petit and Pierre Louis Dulong, _Recherches sur quelques_ _points importants de la th__éorie de la chaleur_, 1819.] It calculates that _any degree of freedom_ (DoF), such as a component $T B D$ of particle translational velocity $T B D$ , which appears _quadratically in energy_, _i.e._ as $T B D$ , contributes $T B D$ to $T B D$ .
For monatomic gases, _e.g._ argon Ar, particle motion is translational only (until atoms dissociate into subatomic particles). With 3 translational DoFs, $T B D$ is constant.

A gas which obeys the ideal gas law with constant $T B D$ is known as _calorically perfect_.

For diatomic gases, _e.g._ nitrogen $T B D$ , molecular motion is translational (left) and rotational (centre) at lower temperature (in the gas state). At higher temperature, vibrational motion is excited along its bond (right). Heat capacity $T B D$ for diatomic gases are then a function of temperature as shown below for $T B D$ .

Rotational motion provides 2 DoFs (rotational energy about the axis is negligible), giving $T B D$ . Beyond 500K, vibrational motion provides 2 additional DoFs (one kinetic, one potential), causing a gradual transition to $T B D$ . At high temperature, molecules dissociate and $T B D$ increases further.

A _thermally perfect_ gas is an ideal gas with $T B D$ . Often $T B D$ can be treated as constant over a range of $T B D$ , _e.g._ T  $T B D$  600K for $T B D$ . Otherwise, accurate calculations require a suitable model of $T B D$ . An _imperfect_ gas exhibits $T B D$ .



== Energy and temperature

Specific internal energy $T B D$ and temperature $T B D$ were described in Sec. 2.17 and Sec. 2.16, respectively. They are related through the specific heat capacity $T B D$ , defined by Eq. (2.61) in Sec. 2.18.

Analyses involving heat usually incorporate both $T B D$ and $T B D$ since:

- $T B D$ is the _measurable_ quantity specified as initial and boundary data and whose data is required as part of the “results”;
- $T B D$ is the _calculated_ quantity solved in energy conservation, _e.g._ Eq. (2.51), but whose data is usually of no interest.


Conversion of values between $T B D$ and $T B D$ is therefore needed, and vice versa. Incorporating Eq. (2.61) into a definite integral for $T B D$ , $T B D$ , gives

$ T B D $
The terms in Eq. (2.62) are illustrated below on a $T B D$ graph. Energy $T B D$ is represented by the area under the curve, in which $T B D$ represents a reference energy up to a reference temperature $T B D$ , and the integral from $T B D$ to $T B D$ is shown by the shaded area.
For applications that cover a reasonably narrow temperature range, $T B D$ can be assumed constant. From Eq. (2.62), the $T B D$ relation becomes

$ T B D $
Alternatively, $T B D$ can be integrated analytically by representing $T B D$ by a polynomial of order $T B D$ with coefficients $T B D$ fitted to measured $T B D$ data
$ T B D $
The values $T B D$ and $T B D$ ultimately add a constant component to $T B D$ . Since Eq. (2.51) is concerned with _changes_ in $T B D$ and the absolute values $T B D$ are usually of no interest, the values of $T B D$ and $T B D$ are often immaterial.
The $T B D$ and $T B D$ values become important when the composition of a fluid changes due to the mixing of constituent fluid species, _e.g._ $T B D$ , $T B D$ , or chemical reactions, _e.g._ with $T B D$ . Each fluid specie possesses a different $T B D$ so any change to the specie concentrations will change $T B D$ of the overall fluid.

In those circumstances, $T B D$ is commonly represented by the _heat of formation_ per unit mass, $T B D$ . The _standard heat of formation_ $T B D$ is the change of enthalpy during the formation of 1 mole of a substance from its constituent elements at standard temperature $T B D$ . Measured heats of formation are available for numerous fluid species.#footnote[Alexander Burcat and Branko Ruscic, _Third Millennium ideal gas and condensed phase thermochemical_ _database for combustion_, 2005.]

If an analysis involves changes to fluid composition, it can then adopt $T B D$ and $T B D$ for individual fluid species, to account for the change in $T B D$ due to changes in the concentrations of fluid species.



== Natural convection

In Sec. 2.13, a set of equations — Eq. (2.47) and Eq. (2.48) — was derived for flow of an incompressible fluid. They are a sample set of equations for mass and momentum conservation that can be solved using methods described in this book.

The example set of equations can be extended to include energy conservation and associated models of heat conduction and heat capacity, described in Sec. 2.15 - 2.18.

The set of equations for mass, momentum and energy can be combined to simulate _natural_ _convection_, _e.g._ for flow of air around a room. In natural convection, a non-uniform temperature causes density variations which generate associated forces due to gravity. Colder air is driven downwards and hot air rises, creating _buoyancy_. Small temperature variations, _e.g._ due to a heat source $T B D$ , can cause buoyancy to be the dominant force.

A simple, approximate equation for $T B D$ can be derived, starting from internal energy conservation in the form of Eq. (2.57). The approximations of constant $T B D$ (with $T B D$ ) and zero viscosity reduce the stress/pressure work terms to zero.

Assuming $T B D$  constant, we can apply Eq. (2.63), which reduces to substituting $T B D$ by $T B D$ since derivatives of constants $T B D$ and $T B D$ are zero. Applying Fourier’s law Eq. (2.54) leads to

$ T B D $
where thermal diffusivity $T B D$ and $T B D$ becomes a thermal source in SI units of $T B D$ .
This is another example, similar to Eq. (2.49), of a transport equation containing a time derivative, advection, diffusion and a source of heat. Applying suitable boundary conditions, the equation can be solved for $T B D$ .

=== Buoyancy force

The effect of buoyancy can be simulated by setting the body force $T B D$ in Eq. (2.47) for an incompressible, Newtonian fluid. While the assumption $T B D$  constant is applied across all the governing equations generally, it cannot be applied to this force. Therefore we apply

$ T B D $
where $T B D$ is a density at a reference state, _e.g._ at the initial $T B D$ and $T B D$ , and $T B D$ is the acceleration due to gravity.
The density $T B D$ is a function of $T B D$ and, optionally, $T B D$ provided by some equation of state, _e.g._ the ideal gas Eq. (2.55). The final momentum equation including this buoyancy force, and assuming $T B D$  constant is:

$ T B D $
The combined set of equations for mass, momentum and energy becomes Eq. (2.48) Eq. (2.67) and Eq. (2.65) respectively, which can be used to solve flows with natural convection.


== Scale similarity

Scale similarity is the notion that, for two systems that are geometrically similar, the flow will follow the same path if the ratio of magnitude of forces acting on the fluid is the same at different points in the flow.

Flows at different scales can be compared using dimensionless variables. The momentum equation, Eq. (2.67), with advection, diffusion and gravity forces can be expressed in non-dimensionalised form by

$ T B D $
with the dimensionless numbers#footnote[After Vincenc Strouhal (1850-1922), Osborne Reynolds (1842-1912), Leonhard Euler (1707-1783) and William Froude (1810-1879).]:
- Strouhal number $T B D$ — transient/steady inertia;
- Reynolds number $T B D$ — inertia/viscous force;
- Euler number $T B D$ — pressure force/inertia;
- (Froude number) $T B D$ $T B D$ — inertia/gravity force.


These dimensionless numbers include a characteristic length $T B D$ , time $T B D$ , speed $T B D$ and pressure $T B D$ . The $T B D$ (hat) notation indicates dimensionless length, time, _etc._, _e.g._ $T B D$ and $T B D$ and corresponding dimensionless operators $T B D$ and $T B D$ .

Equation 2.68 assumes constant $T B D$ and splits gravitational acceleration $T B D$ into its magnitude $T B D$ and unit direction $T B D$ ; pressure, including $T B D$ , is in kinematic units (divided by $T B D$ ).

The dimensionless numbers provide a comparison of the magnitudes of different fluid forces. For example, $T B D$ represents the ratio of inertia force to viscous force and plays a pivotal role in turbulence modelling, introduced in Chapter 6.

Scale similarity applies also to other transported properties. For example, the energy equation, Eq. (2.65), ignoring heat sources $T B D$ , can be expressed in non-dimensional form as

$ T B D $
which includes the following additional dimensionless number:#footnote[After Jean Claude Eugène Péclet (1793-1857).]
- Péclet number $T B D$ — advection/diffusion of heat.


Again in Eq. (2.69), the $T B D$ (hat) notation is applied to temperature $T B D$ to indicate a dimensionless temperature, although it notably does not appear in a dimensionless number.

In fact, with the exception of momentum (which uses $T B D$ ), $T B D$ represents more generally the rates of advection and diffusion, as a ratio, for any transported quantity (here, it is heat).

Further numbers#footnote[After Ludwig Prandtl (1875-1953) and Ernst Heinrich Wilhelm Schmidt (1892-1975).] define the ratios of diffusivities, _e.g._:

- Prandtl number $T B D$ — viscosity/thermal diffusivity;
- Schmidt number $T B D$ — viscosity/mass diffusivity.


where $T B D$ is mass diffusivity (not discussed in this book).

Dimensionless numbers can be multiplied and divided with one another to form further dimensionless numbers. For example, the Péclet number for heat transfer $T B D$ .



== Region of influence

A CFD calculation is performed by solving partial differential equations, such _e.g._ for conservation of mass, momentum and energy, over a solution domain. It requires suitable boundary conditions, discussed in Chapter 4.

The solution is influenced by any change in a field value, _e.g._ $T B D$ , $T B D$ , $T B D$ , at some point in the domain, _e.g._ at an inlet boundary. The form of the equation determines the way in which these changes, or _perturbations_, propagate across the domain over time.

Momentum and mass conservation for an incompressible fluid can be represented by Eq. (2.67) and Eq. (2.48). The form of Eq. (2.48), including only a Laplacian derivative $T B D$ , ensures that $T B D$ is influenced instantaneously at all points in the domain by a perturbation in $T B D$ at any point.

The resulting instantaneous change in $T B D$ then causes $T B D$ to be redistributed everywhere by Eq. (2.67), with further short-range changes due to advection and diffusion (discussed next).

The outcome for an incompressible fluid is that $T B D$ is influenced instantaneously everywhere in the domain by a perturbation at any point. In other words, the speed of sound $T B D$ , corresponding to propagation of disturbances, is infinite.

=== Advection-diffusion equations

Energy conservation can be represented by Eq. (2.65), in which perturbations propagate at a characteristic speed due to the advection and diffusion.

Advection propagates at speed $T B D$ , with a region of influence $T B D$ in the direction of flow. The relation can be obtained from scale similarity when $T B D$ and $T B D$ are similar in magnitude, _i.e._ $T B D$ in Eq. (2.69).

By the same argument, the region of influence for diffusion coincides with $T B D$ and $T B D$ being similar in magnitude, _i.e._ $T B D$ . A diffusion “front” travels a distance according to $T B D$ ( $T B D$ means “of the order of magnitude”).

The diffusion equation, $T B D$ in one dimension ( $T B D$ ) has a solution#footnote[Horatio Carslaw and John Jaeger, _Conduction of Heat in Solids_.] $T B D$ where $T B D$ , $T B D$ are constants and $T B D$ is the error function below.

This means that for heat conduction problems, _e.g._ in solids, the distance travelled by the thermal front is $T B D$ , consistent with similarity arguments above.

The coefficient $T B D$ depends on where on the $T B D$ curve the front is located. One option is $T B D$ , where the solution is within 0.5% of the asymptote of 1.



== Summary of equations

- *Force* $T B D$ at a surface:
$ T B D $
- *Rate of deformation*
$ T B D $


=== Conservation laws

- *Mass*, Sec. 2.4
$ T B D $
- *Momentum*, Sec. 2.7
$ T B D $
- *Energy*, Sec. 2.15
$ T B D $


=== Constitutive models

- *Newtonian fluid*, Sec. 2.12
$ T B D $
- *Fourier’s law*, Sec. 2.16
$ T B D $
- *Ideal gas*, Sec. 2.16
$ T B D $
- *Specific heat capacity*, Sec. 2.18
$ T B D $
- *Energy-temperature relation*, Sec. 2.19
$ T B D $


=== Derivatives

- *Material time derivative*, Sec. 2.5
$ T B D $
- *Divergence*, Sec. 2.4, $T B D$ is vector, tensor
$ T B D $
- *Gradient*, Sec. 2.12
$ T B D $
- *Laplacian*, Sec. 2.14, $T B D$
$ T B D $




== Summary of tensor algebra

Below is tensor algebra applied to Cartesian ( $T B D$ ,  $T B D$ ,  $T B D$ ) co-ordinates using: scalar $T B D$ ; vectors $T B D$ , $T B D$ ; tensors $T B D$ , $T B D$ , $T B D$ .

=== Products

- *Inner product of two vectors*, Sec. 2.3
$ T B D $
- *Outer product of two vectors*, Sec. 2.8
$ T B D $
- *Inner product of vector and tensor*, Sec. 2.6
$ T B D $
- *Inner product of two tensors*, Sec. 2.8
$ T B D $
- *Double inner product of two tensors*, Sec. 2.17
$ T B D $
- *Cross product of two vectors*, first used in Sec. 3.3, produces a vector with components
$ T B D $


=== Tensors and operations

- *Transpose of a tensor*, Sec. 2.7
$ T B D $
- *Symmetric and skew tensors*, Sec. 2.7
$ T B D $
- *Trace of a tensor*, Sec. 2.10
$ T B D $
- *Identity tensor*, Sec. 2.8
$ T B D $
- *Deviatoric and spherical tensors*, Sec. 2.10
$ T B D $




== Vector identities

Several identities are listed below which can be verified under the assumption that the relevant derivatives exist and are continuous. The identities are expressed for: scalars $T B D$ , $T B D$ ; vectors $T B D$ , $T B D$ , $T B D$ ; tensor $T B D$ .

=== Algebraic operations

$ T B D $

=== Gradient derivatives

$ T B D $

=== Divergence derivatives

$ T B D $

=== Laplacian derivatives

$ T B D $

=== Curl derivatives

$ T B D $

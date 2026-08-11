Chapter 3 described:

- the computational mesh consisting of cells of any irregular polyhedral shape with mesh properties for the finite volume method;
- discretisation of terms in equations to construct a matrix equation $T B D$ for solution variable $T B D$ ;
- discretisation schemes for optimal accuracy, boundedness and convergence, using limiting and correction;
- the time derivative, Courant number and time step size.


At the boundary of the mesh, constraints must be applied to the solution variable $T B D$ , known as _boundary conditions_. Setting boundary conditions is challenging because:

- they need to _reflect the physical conditions_ of the case being simulated;
- they need to make the case _well posed_, _i.e._ provide unique, stable, physical solutions to the equations;
- they must be compatible across sets of multiple equations, in particular $T B D$ and $T B D$ .


The mesh boundary is split into regions known as _patches_, on which different boundary conditions are applied. The choice of boundary condition generally depends on flow direction at a patch, whether the patch corresponds to a solid wall, _etc._

There are basic forms of boundary condition which specify the value, normal gradient, _etc._ of $T B D$ at the boundary. They are applied through modifications to the matrix coefficients $T B D$ and source vector $T B D$ using mesh data of the faces and cells adjacent to the boundary.

The opening topics of this chapter describe the mesh data and the numerical methods required by the basic forms of boundary condition.

More specialised boundary conditions can be _derived_ from the underlying forms which correspond to different physical conditions. Some derived conditions, that are often used for particular boundary configurations, are introduced in this chapter.

Otherwise, the number of possible derived conditions is almost unlimited due to the range of potential physical conditions that can be encountered in fluid flow problems. It is left to the producers of CFD software to document those conditions, when it is important that the user identifies the basic underlying type of the condition, _e.g._ fixed value or gradient.



== Boundary mesh

Sec. 3.2 describes the computational mesh used by the finite volume method. It identifies a domain boundary by faces which are connected only to one cell. The boundary faces are grouped into _patches_, each under a unique name, corresponding to different regions of the boundary. Boundary conditions are then applied to each patch to provide a representative solution to the case of interest.

In the example of flow in a pipe, it would be logical to split the boundary into 3 patches in order to specify inflow, outflow and solid walls; and, _inlet_, _outlet_ and _wall_ would be logical names for these patches.

=== Patch geometric data

The geometry of a patch is described using face data described in Sec. 3.3, including:

- the face area vector $T B D$ , with area magnitude and direction  $T B D$ ;
- the face unit normal vector $T B D$ ;
- the face centre $T B D$ .


The cell connected to each patch face is denoted by the subscript “P”, _e.g._ the cell centre is denoted $T B D$ .

=== Patch deltas

Sec. 3.8 describes the “delta” $T B D$ for each face as the vector connecting the centres of its owner and neighbour cells. $T B D$ is defined differently for a patch face since it has no neighbour cell.

The delta is defined as the component of $T B D$ in the face normal direction. The surface gradient $T B D$ is orthogonal to the face, which eliminates the error associated with non-orthogonality, at the expense of introducing a skewness error. Taking the inner product with $T B D$ gives the magnitude which is then multiplied by $T B D$ to assign the direction, _i.e._

$ T B D $
The _delta coefficients_ are $T B D$ , as defined in Sec. 3.8. The delta coefficient $T B D$ , representing “inverse distance”, is a critical parameter in the discretisation of boundary conditions.


== Fixed value and fixed gradient

In an equation for $T B D$ that is discretised to form a matrix equation $T B D$ , there are two terms that include properties interpolated to faces:

- an advection term of the form $T B D$ which is discretised by Eq. (3.8) in extensive form (see Sec. 3.5) as $T B D$ ;
- a Laplacian term of the form $T B D$ which is discretised by Eq. (3.2) in extensive form as $T B D$ .


The advection term requires the value $T B D$ and Laplacian term requires the surface normal gradient $T B D$ at faces. When a face is part of a boundary patch, any gradient $T B D$ and/or value $T B D$ must be specified through boundary conditions.

The *fixed value*, or _Dirichlet_ condition,#footnote[Peter Dirichlet, _Über einen neuen Ausdruck zur Bestimmung der Dictigkeit einerunendlich d__ünnen_ _Kugelschale, wenn der Werth des Potentials derselben in jedem Punkte ihrer Oberfl__äche gegeben ist_, 1850.] is the first type of boundary condition, where the boundary value $T B D$ is specified. For example, at an inlet patch, we might specify a temperature $T B D$  K of the fluid flowing into the domain.

The *fixed gradient*, or _Neumann_ condition,#footnote[after Carl Neumann 1832-1925.] is the second type, in which the gradient normal to the boundary $T B D$ is specified (where $T B D$ ). In many cases, the applied normal gradient is zero, which is a common condition applied to many fields, including $T B D$ , at an outlet patch, _i.e._ where the fluid flows out of the domain.

When an advection term is discretised, a fixed value condition is applied by substituting the face value $T B D$ with the patch value $T B D$ , _i.e._ setting $T B D$ . When a fixed gradient $T B D$ is specified, the face value is expressed as follows, where $T B D$ is the value in the cell adjacent to each face:

$ T B D $
When a Laplacian term is discretised, a fixed gradient boundary condition is applied by substituting the face normal gradient $T B D$ with the patch face normal gradient $T B D$ , _i.e._ setting $T B D$ . When a fixed value $T B D$ is specified, the face normal gradient is expressed by
$ T B D $



== Fundamentals of boundary conditions

The specification of boundary conditions is one of the most challenging tasks in setting up a CFD simulation. The range of possible boundary conditions is endless, to cover all of the potential applications and physics.

Setting boundary conditions is not an exact science but is guided by an basic specification using the fixed value (Dirichlet) and fixed gradient (Neumann) conditions introduced in Sec. 4.2.

The figure above shows the basic specification for $T B D$ , $T B D$ and $T B D$ for incompressible subsonic flow, _e.g._ described by Eq. (2.47), Eq. (2.48) and Eq. (2.65).

The boundary conditions for $T B D$ and $T B D$ require particular attention since they are coupled. The conditions on $T B D$ are independent of $T B D$ and $T B D$ and are representative of other transported scalar fields, _e.g._ turbulent kinetic energy $T B D$ .

The specification for inlets and outlets summarises as:

- zero gradient on $T B D$ at an inlet, fixed value on other variables;
- fixed value on $T B D$ at an outlet, zero gradient on other variables.


=== Propagation of disturbances

A disturbance in a flow is simply any change from an equilibrium or steady solution. A disturbance at one location travels, or _propagates_, through the fluid.

The combination of boundary conditions at open boundaries, _i.e._ those excluding walls, relates to the propagation of disturbances. While disturbances are transported by advection with the flow, they propagate as waves at the speed of sound $T B D$ .

Sound waves can propagate disturbances against the direction of flow if it is _subsonic_, _i.e._ $T B D$ . Disturbances must be able to propagate _outwards_ through the inlet, which requires one variable not to be prescribed at the inlet. Similarly, they must be able propagate _inwards_ from the outlet, which requires one variable to be prescribed at the outlet.

A pressure equation, combining mass and momentum conservation, describes wave propagation. For an incompressible fluid the wave speed $T B D$ , since Eq. (2.48) contains no $T B D$ term. A disturbance at any point influences the solution everywhere in the domain instantaneously, as discussed in Sec. 2.22.

Pressure is then logically the variable on which to specify the boundary conditions to support wave propagation. Pressure is therefore: prescribed at the outlet, _i.e._ we specify fixed value condition; not prescribed at the inlet, _i.e._ we specify a fixed gradient condition, usually set to $T B D$ .



== Wall boundaries

The *no-slip* condition is generally applied at a solid wall which is impermeable (assuming $T B D$ ). The condition is $T B D$ , where $T B D$ is the velocity of the wall, which is usually stationary with $T B D$ . The proof behind the no-slip condition is that it predicts a pressure drop along tubes of small diameter which matches experiments.#footnote[Jean Poiseuille, _Recherches exp__érimentales sur le mouvement des liquides dans les tubes de tr__ès petits_ _diam__ètres I-IV_, 1840.]

The 2D, lid-driven cavity is a flow problem in which no-slip conditions are applied at all boundaries. It provides insight into the boundary condition for $T B D$ at a wall. From Eq. (2.47) for an incompressible fluid, with $T B D$  constant, the component normal to the domain boundary is

$ T B D $
At the solid wall boundary, the last 2 terms in Eq. (4.4) disappear since $T B D$ , reducing the gradient condition to
$ T B D $
The second term can be written ( $T B D$  constant) $T B D$ , so is only non-zero where there is flow normal to the boundary in its vicinity, _e.g._ at the corners of the cavity in our example. The term is usually small and its calculation involves extrapolation from the internal solution which often causes instability, so it is generally ignored.
A body force $T B D$ , _e.g._ gravity, is generally prescribed so it does not introduce instability. Where it is significant, it must be included in the boundary condition, _i.e._ $T B D$ .

Otherwise, in the absence of a body force, we reach the standard form of boundary condition for pressure at a wall, $T B D$ .

=== Fixing pressure

With only fixed _gradient_ conditions on pressure at the boundary, the pressure _value_ is not fixed at any point in solution domain. The solution is not unique, as shown in the 1D example below with gradient conditions at both ends.

To achieve a unique solution, $T B D$ must then be fixed to a reference value $T B D$ at a reference cell $T B D$ in the domain. To achieve this, the diagonal coefficient $T B D$ is doubled and $T B D$ is added to source $T B D$ , in the matrix equation $T B D$ described by Eq. (3.1). This minimal change “pins” the solution to $T B D$ in cell $T B D$ .



== Inlets and outlets

A basic set of boundary conditions was introduced in Sec. 4.3 for incompressible subsonic flow.

At an _inlet_, fields are generally specified as fixed value. From a “physical” perspective, this is justified by the fact that disturbances propagate in the direction of flow so must be specified at the upstream boundary.

Advection requires interpolation of variables from cell centres ( $T B D$ ) to faces ( $T B D$ ). Some degree of upwind interpolation, _e.g._ as part of a limited scheme, is generally required. At a face at an inlet patch there is no upwind cell, so values $T B D$ must be specified at inlet faces instead. From a “numerical” perspective, this justifies the fixed value boundary condition at an inlet.

For incompressible, subsonic flow, $T B D$ is the exception. A fixed gradient condition allows disturbances to propagate upstream through the inlet as sound waves. The gradient condition is further justified since there is no advection of $T B D$ , see Eq. (2.48), so no upwind interpolation is required at the inlet patch faces.

At an _outlet_, the converse is then true: fields are specified as fixed gradient condition, with the exception of $T B D$ which is fixed value. The outlet conditions ultimately dictate the traction force combining Eq. (2.16), Eq. (2.33), Eq. (2.46) and Eq. (2.41) as follows:

$ T B D $

The standard condition applied to $T B D$ is $T B D$ . By Eq. (4.6), this results in a uniform _normal_ traction force corresponding to the outlet pressure $T B D$ . The traction force _tangential_ to the outlet $T B D$ , where $T B D$ is the tangential direction and $T B D$ is the normal component of $T B D$ .

=== Supersonic conditions

If the fluid is compressible and the flow speed is _supersonic_ at the inlet, _i.e._ $T B D$ , waves can no longer propagate outwards through the inlet boundary. When this occurs, a fixed value $T B D$ must be specified at the inlet.

Similarly if the flow is supersonic at an outlet, all disturbances propagate through the outlet. In that case, $T B D$ cannot be specified, _i.e._ a condition on the normal gradient $T B D$ is applied.

Some flow domains combine subsonic flow at the inlet and supersonic flow at the outlet. When this occurs, a gradient condition is required for $T B D$ at both boundaries, leaving $T B D$ under-specified. This problem is best overcome by moving the outlet boundary sufficiently downstream for the flow to expand to subsonic speed. This allows $T B D$ to be specified through a fixed value condition.



== Free (entrainment) boundaries

The basic set of boundary conditions for subsonic incompressible flow, introduced in Sec. 4.3, relates to inlet, outlet and wall boundary regions. Sometimes a single boundary patch occupies a location which experiences both inflow and outflow.

The figure above shows a 2D domain with an inlet and outlet, left and right respectively. The top boundary is _free_ to allow:

- inflow from the centre to the left side of the patch, through entrainment of fluid by the main flow from the inlet;
- outflow at the right side of the patch, driven by a pressure gradient that emanates from the point where the flow impinges on the right wall.


Boundary conditions for $T B D$ and $T B D$ are therefore required which support both inflow and outflow at a single patch. The inlet flow speed is determined from the flow solution rather than being prescribed at the patch. This suggests a gradient condition on $T B D$ , which is generally prescribed at an outlet.

At a free/entrainment boundary, the basic outflow conditions, $T B D$ and $T B D$ , are generally *not* recommended. The figure above shows a snapshot from a simulation using those conditions. The solution oscillates between changing levels of inflow and outflow over the open boundary; the figure below shows the flow leading up to the snapshot above.

The is no convergence towards a stable solution in this example of a simple flow in a regular 2D geometry. Moreover, solutions tend to diverge and “blow up” if the basic outflow conditions are applied at a free boundary for geometries and meshes that are even moderately complex. In the following sections we look at some alternative conditions that are applied at a free boundary to deliver stable, robust solutions.



== Total pressure condition

Sec. 4.6 concluded that the basic outflow conditions — $T B D$ and $T B D$ — are generally unstable for a free boundary with both inflow and outflow.

The conditions are unstable due to pressure fluctuations at the boundary, shown above. Flows oscillates in and out, shown by 3 stages, creating a vortex that travels from left to right:

+ the pressure gradient, decreasing outward, causes outflow;
+ the outflow speed increases causing the pressure gradient to change direction;
+ inflow begins and the speed increases, until the pressure gradient changes direction, returning back to step 1.


The *total pressure* boundary condition improves the stability of solutions. It is a fixed value type, calculated according to:

$ T B D $
The specified total pressure, $T B D$ , can be imagined as the fluid pressure under quiescent conditions far from the free boundary, which decreases as the fluid accelerates towards the boundary. Note that Eq. (4.7) is written for the incompressible assumption, Eq. (2.47), where $T B D$ and $T B D$ are kinematic, _i.e._ divided by $T B D$ .
The solution using the total pressure condition converges to the flow field shown above. The critical effect of this boundary condition is that, the boundary $T B D$ decreases by $T B D$ as the inflow speed $T B D$ increases. This reduces the pressure gradient driving inflow, which moderates the increase in inflow speed, enabling it to settle to a stable level.

=== Total pressure for high speed flow

The total pressure condition can be applied to high-speed flow of a compressible gas. The calculation of $T B D$ for inflow is simply replaced with the 1D isentropic flow equation,#footnote[Ascher Shapiro, _The dynamics and thermodynamics of compressible fluid flow, Vol. 1, Ch. 4_, 1953.]

$ T B D $
where Mach number is $T B D$ and $T B D$ .


== Numerical framework

The fixed value and fixed gradient conditions described in Sec. 4.2 can be combined to form a general numerical framework for boundary conditions.



The contributions from the boundary conditions to the matrix equation $T B D$ , by discretisation of advection and Laplacian terms, can be generalised as:

- “_internal_” contributions to $T B D$ , from terms including the cell value $T B D$ ;
- “_boundary_” contributions to $T B D$ , from terms without $T B D$ .


The example above shows Eq. (4.2) for the face value, required at a boundary for advection discretisation, in the case of a _fixed gradient_ boundary condition. The internal “factor” on $T B D$ is 1, which is multiplied by $T B D$ for the contribution to the respective diagonal coefficient in $T B D$ , as in the example in Sec. 3.24.

The boundary factor is $T B D$ , which is similarly multiplied by $T B D$ for the contribution to $T B D$ .

For the _fixed value_ condition $T B D$ with _advection_, the boundary factor is $T B D$ and an internal factor is 0.

Laplacian discretisation requires the surface normal gradient on the faces. A _fixed gradient_ condition delivers an equivalent boundary factor of $T B D$ to $T B D$ and an internal factor of 0.


With a _fixed value_ $T B D$ , the face normal gradient Eq. (4.3) gives an internal factor of $T B D$ and a boundary factor of $T B D$ . Both are multiplied by $T B D$ in their contributions to $T B D$ and the diagonal coefficient of $T B D$ in $T B D$ , as shown in the Laplacian discretisation in Sec. 3.24.
The table below summarises: the “value” internal and boundary factors, contributing to the respective matrix coefficients with advection discretisation; and equivalent “gradient” factors relating to Laplacian discretisation. This provides a framework which can be extended to more complex conditions.


term factor fixed value fixed gradient advection value internal 0 1 advection value boundary $T B D$ $T B D$ Laplacian gradient internal $T B D$ 0 Laplacian gradient boundary $T B D$ $T B D$


== Mixed fixed value/gradient

Sec. 4.8 concluded with a table of factors for the contributions to coefficients of $T B D$ and $T B D$ for the fixed value and fixed gradient boundary conditions.

It distinguishes between contributions for the discretisation of an advection term which requires _values_ at faces, and the Laplacian term which requires the normal _gradient_.

A *mixed* fixed value/gradient condition is defined by introducing a _value fraction_ $T B D$ for which

$ T B D $
Within the range $T B D$ the condition operates in between fixed value and fixed gradient. The mixed condition is simply implemented by blending the fixed value and fixed gradient contributions by $T B D$ , as shown in the table below.

factor mixed value internal $T B D$ value boundary $T B D$ gradient internal $T B D$ gradient boundary $T B D$
This mixed condition provides the framework for a boundary condition that can switch between a fixed value $T B D$ and fixed gradient $T B D$ , by changing $T B D$ . Switching is often based on flow direction, $T B D$ corresponding to inflow and $T B D$ to outflow.

Some boundary conditions can operate in the range of value fractions $T B D$ . The Robin condition, described next, can also be expressed as a mixed condition with varying $T B D$ .

=== Robin condition

The _Robin_ condition#footnote[after Victor Robin (1855-1897).] combines the value and normal gradient at the boundary through an expression:

$ T B D $
where $T B D$ is a scalar coefficient with units of length and $T B D$ is some constant value of $T B D$ .
The Robin condition can be treated like the mixed condition by relating $T B D$ to a reference fixed value $T B D$ and gradient $T B D$ at a boundary, according to $T B D$ .

Substituting for $T B D$ in Eq. (4.10) and making $T B D$ the subject of the equation gives:

$ T B D $
Comparison with the value factors in the previous table shows that the Robin condition corresponds to the mixed condition with the value fraction $T B D$ .
In this form, $T B D$ and $T B D$ relate to the limits $T B D$ and $T B D$ , respectively. Values can be selected in these limits to represent the physics of the condition.

In many cases the reference gradient is $T B D$ such that $T B D$ in Eq. (4.10). For example a condition for temperature $T B D$ would tend to $T B D$ as $T B D$ and $T B D$ as $T B D$ .

The value fraction $T B D$ includes $T B D$ , so the condition operates “in the middle” between the fixed value and gradient when $T B D$ is the same order of magnitude as the boundary cell height.



== Mixed inlet-outlet condition

The *inlet-outlet* boundary condition is the most basic example of the mixed fixed value/gradient type, described in Sec. 4.9. The condition sets the reference gradient $T B D$ and uses a specified reference value $T B D$ . The value fraction is then set to

$ T B D $
The flow direction is established from the sign of the volumetric or mass flux $T B D$ at each boundary face, described in Sec. 2.8, by
$ T B D $
The inlet-outlet condition is generally very useful for _scalar_ fields, _e.g._ turbulence fields, $T B D$ , _etc._ It has an immediate practical use at a free boundary, _e.g._ in the case introduced in Sec. 4.6.
The figure shows the solution of Eq. (2.65), converged over time with $T B D$ and unity Prandtl number $T B D$ , see Sec. 2.21. The fixed condition $T B D$ is applied at the inlet and a zero gradient condition $T B D$ at the walls.

At the free boundary, the inlet-outlet condition enables $T B D$ to be specified where inflow occurs. The inlet value in the example is set to $T B D$ ; the image shows mixing of fluids at different temperatures, from the inlet and entrained at the free boundary.

=== Numerical benefit of inlet-outlet

Boundaries may be described “inlet” and “outlet” based on the expectation of the flow direction during a simulation. But the flow direction may not always happen as expected.

In the case of an outlet, for example, inflow might occur during a simulation. For example, at the start of a simulation, the initial conditions may induce inflow before the internal flow is established. Localised inflow can also occur when rotating flow structures pass through an outlet boundary, _e.g._ when a bluff body sheds vortices, as shown below.

Where inflow occurs, the inlet-outlet condition can switch to the fixed value type to ensure stability, as discussed in Sec. 4.5. The inlet-outlet condition is therefore commonly applied to scalar fields (except $T B D$ ), at a boundary which is notionally an outlet, to avoid numerical instability associated with unexpected inflow.



== Transform condition

Some boundary conditions represent $T B D$ at the boundary as a transformation of the cell value $T B D$ . They can be expressed in terms of a general *transform* condition

$ T B D $
where $T B D$ is the geometric transformation of variable $T B D$ by a tensor $T B D$ . The transformation is calculated as follows:
$ T B D $
When $T B D$ is a scalar, the transform condition is equivalent to a zero gradient condition. Otherwise, it is implemented so that terms in $T B D$ in Eq. (4.13) contribute to coefficients in $T B D$ . These contributions are implicit which improves convergence when solving the matrix equation $T B D$ .
A factor $T B D$ is introduced to specify the contribution to internal coefficients. It represents a single internal coefficient for each component of $T B D$ so has the same rank as $T B D$ , _i.e._ it is a vector when $T B D$ is a vector, and a tensor when $T B D$ is a tensor.

The multiplication of each coefficient of $T B D$ by its respective component of $T B D$ is denoted by by $T B D$ . In the case of vector $T B D$ , this “component multiplication” is

$ T B D $
For advection discretisation, the face value is represented as
$ T B D $
where $T B D$ is an explicit boundary value, _calculated_ from the expression in Eq. (4.13) using the current values $T B D$ . The other explicit term uses the current $T B D$ with “ $T B D$ ” denoting $T B D$ for a vector $T B D$ and $T B D$ for a tensor.
Laplacian discretisation requires the face normal gradient $T B D$ . Combining $T B D$ with Eq. (4.15) gives

$ T B D $
where the explicit gradient is calculated by
$ T B D $
The transform condition is summarised within the table below by the value and gradient contributions.

factor transform value internal $T B D$ value boundary $T B D$ gradient internal $T B D$ gradient boundary $T B D$


== Symmetry condition

The transform boundary condition was presented in Sec. 4.11. It provides a convenient framework for implementing boundary conditions that represent a _geometric constraint_, including the _symmetry_ condition.

The *symmetry* condition is suitable for simulations where the geometry contains a plane of symmetry and the flow field is assumed symmetric. By generating a mesh on one side of a plane of symmetry and applying the symmetry condition, the number of cells, and hence solution time, is reduced.

In the context of a wall boundary, the symmetry condition is also equivalent to *slip* (as opposed to the common no-slip condition).

A symmetry plane is a transform condition so when the solution variable $T B D$ is a scalar, it reduces to zero gradient. For a vector, _e.g._ $T B D$ , the condition is zero gradient tangential to the plane, and zero fixed value normal to the plane.

When $T B D$ is a tensor, the boundary condition requires a more precise definition, which can also be applied to a vector $T B D$ . The boundary $T B D$ can be considered as the mean of the adjacent cell $T B D$ and the mirror image $T B D$ transformed by the reflective transformation tensor $T B D$ , _i.e._

$ T B D $
Here, $T B D$ is the unit normal vector on the boundary face.
Using the notation in Sec. 4.11, the explicit boundary value $T B D$ is calculated using current $T B D$ from Eq. (4.18) and the explicit gradient $T B D$ is calculated by Eq. (4.17).

Comparing Eq. (4.18) with the transform condition of Sec. 4.11, the factor $T B D$ corresponds to the tensor $T B D$ . For a vector field, a factor which gives good solution convergence is

$ T B D $
where “ $T B D$ ” is the vector of the diagonal components of tensor $T B D$ .
For a tensor field, good convergence is achieved with a tensor $T B D$ calculated as the outer product of $T B D$ for a vector, _i.e._ denoting the vector factor by $T B D$ , the tensor factor is $T B D$ .

=== Orthogonality condition

The axes $T B D$ , $T B D$ , $T B D$ , introduced in Sec. 2.1, must remain orthogonal under a transformation. It requires the transpose of the transformation tensor $T B D$ to equal its inverse, _i.e._ $T B D$ .

The orthogonality condition is therefore $T B D$ The reflective transformation $T B D$ satisfies the orthogonality condition since

$ T B D $



== Axisymmetric (wedge) condition

There are some fluid flow problems for which the geometry is axisymmetric. Assuming the flow solution is axisymmetric, _i.e._ fields do not change in the circumferential direction, the computational mesh can be formed of a wedge-shaped slice of the flow geometry.

This type of mesh for axisymmetric solution contains one cell across the circumferential direction, which reduces the number of cells to two dimensions in the axial and radial directions.

This approach to axisymmetric solution introduces a geometric error due the faces normal to the radial direction being flat. This error reduces with decreasing wedge angle; in practice, the error can be considered negligible for an angle of 1 $T B D$ .

The *wedge* boundary condition is applied to the two sloping side patches. It transforms cell values $T B D$ to the patch faces using a rotational transformation tensor $T B D$ by

$ T B D $
$T B D$ defines a rotation between the unit vector $T B D$ in the circumferential direction at the cell centre and the unit face normal vector $T B D$ by $T B D$
The wedge condition uses the general transform framework of Sec. 4.11, with the explicit value $T B D$ calculated using current $T B D$ from Eq. (4.20). The explicit gradient $T B D$ is the boundary gradient $T B D$ calculated from $T B D$ in an imaginary neighbour cell by

$ T B D $
where the rotation between cell centres $T B D$ .
The factor $T B D$ is chosen to minimise the gradient boundary coefficients (see Sec. 4.11). The vector factor is $T B D$ where “ $T B D$ ” is defined in Sec. 4.12.

For a tensor $T B D$ , the factor is $T B D$ , where $T B D$ are the coefficients for a vector, as described for the symmetry condition in Sec. 4.12.

=== Rotation tensor

The rotation tensor $T B D$ between two _unit_ vectors $T B D$ and $T B D$ can be calculated using the Euler-Rodrigues rotation formula,#footnote[Benjamin Olinde Rodrigues, _De l’attraction des sph__éro__ïdes_, 1815.]

$ T B D $
where $T B D$ and $T B D$ .


== Direction mixed condition

In some cases, a boundary condition is sometimes needed which applies a different underlying type — fixed value or gradient — to different components of a non-scalar field.

The condition is most readily applied to velocity $T B D$ , where different conditions are applied to its normal and tangential components to the patch, $T B D$ and $T B D$ respectively.

The *direction mixed* condition combines the mixed condition from Sec. 4.9 with the transform conditions of Sec. 4.11. The mixed condition can be first expressed in the form of a general transform condition Eq. (4.13) as

$ T B D $
In the direction mixed condition the value fraction $T B D$ is replaced by a transformation tensor $T B D$ , whose components are value fractions in the range 0 to 1, by
$ T B D $
The value fraction tensor $T B D$ is set according to the requirements of the boundary condition which is derived from this direction mixed framework.
Imagine an example of $T B D$ at a face oriented with normal vector $T B D$ , for which the $T B D$ condition is fixed gradient and $T B D$ is fixed value. The value fraction must be 1 in the tangential direction and 0 in the normal direction, which gives:

$ T B D $
The boundary condition is then implemented similarly to the symmetry and wedge conditions of Sec. 4.12 and Sec. 4.13. First, $T B D$ is the calculated $T B D$ using current $T B D$ by Eq. (4.23) and the explicit gradient $T B D$ is calculated from Eq. (4.17).
The factor $T B D$ corresponds to $T B D$ . For a vector field, it is calculated, as in the symmetry condition in Eq. (4.19), by

$ T B D $
For a tensor $T B D$ , the factor is $T B D$ , where $T B D$ are the coeffcients for a vector, as described for the symmetry condition in Sec. 4.12.
The condition is implemented using value and gradient factors according to the transform condition, summarised in the table on page 284. Any boundary condition which is based on this direction mixed condition then only requires a description of the $T B D$ , $T B D$ and $T B D$ parameters.



== Inlet-outlet-velocity condition

The inlet-outlet condition was described in Sec. 4.10, which switches between the zero gradient type for outflow and fixed value for inflow. The condition is not so suitable for velocity, _e.g._ at a free boundary, since it is requires an inlet velocity to be _prescribed_ when inflow occurs.

The switching based on the flow direction is not a problem in itself. The flow direction comes from the sign of the flux $T B D$ , as noted in Eq. (4.12), which is determined by the solution of the pressure equation (part of pressure-velocity coupling described in Chapter 5).

The problem is instead that a value of velocity cannot be chosen in the case of inflow since the inflow speed is _determined_ by the solution within the domain. This suggests a zero gradient type as a more suitable condition.

However, Sec. 4.3 discusses that all but one variable should be prescribed at a boundary with inflow. For velocity that variable could be its _normal_ component, while a fixed value could still be applied to the _tangential_ component.

The direction mixed type described in Sec. 4.14 provides the framework to implement this condition. A reference value $T B D$ is specified and the reference gradient is zero, _i.e._ $T B D$ .

The value fraction tensor is calculated as

$ T B D $
The fixed value which is applied must be a velocity $T B D$ tangential to the boundary, which can be calculated by subtracting the normal component from a specified $T B D$ by $T B D$
This *inlet-outlet-velocity* condition can be applied at the free boundary in the example from Sec. 4.6. For free boundaries like this, $T B D$ is the only practical specification, which causes all inflow to be normal to the free boundary, as shown above.

The solution is clearly different to that using $T B D$ shown on page 267. There, the inflow direction is determined by the solution rather than prescribed. The solution with $T B D$ may be more accurate but the inlet-outlet-velocity condition adheres better to the general principles of boundary conditions, so is likely to be more stable.

It becomes less contentious to set $T B D$ with the inlet-outlet-velocity condition as the boundary is moved further into the far field, where the flow is more quiescent.



== Blended freestream condition

There is a class of problems in CFD that involve external flow around one or more solid bodies, _e.g._ a vehicle, wind turbine, buildings, _etc._

A solution domain is specified which includes the solid body and extends some distance to a free boundary in the _far-field_. A flow velocity $T B D$ is specified which can be applied as a fixed value type at an inlet patch.

The far-field boundary requires attention. The robust conditions at a free boundary for $T B D$ and $T B D$ are inlet-outlet-velocity and total pressure described in Sec. 4.15 and Sec. 4.7 respectively.

The inlet-outlet-velocity requires $T B D$ to be prescribed for inflow which may differ significantly from a determined $T B D$ when there is outflow. Solution accuracy depends on the suitability of the prescribed $T B D$ .

The flow direction in the far-field can often be close to tangential to the boundary, especially with a box-shaped domain. If the flow at one face changes from outflow to inflow, $T B D$ suddenly changes to the prescribed fixed value and $T B D$ decreases by $T B D$ . Sometimes a pattern of switching can occur in adjacent faces and repeated switching can slow the convergence of a solution.

The *blended freestream* condition is a mixed type with zero reference gradient, $T B D$ , which modifies the value fraction $T B D$ as shown above. In the limit that the flow direction is normal to the boundary, the condition becomes the fundamental fixed value and zero gradient types for $T B D$ and $T B D$ .

Between these extremes, $T B D$ is blended linearly, _e.g._ for $T B D$ by

$ T B D $
This means $T B D$ for both $T B D$ and $T B D$ for tangential flow.
At a boundary face, $T B D$ may be directed normal-inward, causing $T B D$ by Eq. (4.27). The condition can then “lock” at $T B D$ , so to avoid this, the calculation can use a velocity equating to the mean of the face and neighbour cell value, _i.e._ $T B D$ .

Note that for $T B D$ , the value fraction is calculated changing the sign of the second term in Eq. (4.27) _i.e._ $T B D$

The freestream conditions overcome the problem of switching to improve the convergence of solutions. Boundary velocities are determined, not prescribed, which seems to improve accuracy, _e.g._ in force calculations described in Sec. 8.4 and Sec. 8.6.



== External wall heat flux

The boundary condition for temperature (or energy) dictates the heat transfer across a boundary. At a boundary that represents a solid wall, simple conditions can sometimes be applied. However, specialised boundary conditions are often required that control the heat flux across the boundary.

The *fixed temperature* is the simplest condition, setting a fixed value $T B D$ . This condition provides an approximation for cases for a solid with high _thermal mass_, due to a large mass of material and high conductivity $T B D$ , which helps to maintain constant $T B D$ .

Otherwise the boundary condition sets the heat flux normal to the boundary, $T B D$ derived from Eq. (2.54) by

$ T B D $
Another simple condition is zero gradient $T B D$ . This is the *adiabatic* condition, corresponding to zero normal heat flux by Eq. (4.28), suitable when the solid is a thermally insulating material with a large mass and low $T B D$ .
Otherwise a *fixed heat flux* condition specifies an _inward_ heat flux $T B D$ as a fixed gradient type with a reference gradient by

$ T B D $

=== Fixed heat transfer coefficient

Another way to specify the heat transfer at an external wall is by Newton’s law of cooling.#footnote[Isaac Newton, _Scala graduum caloris. Calorum descriptiones & signa_, Philosophical Transactions, *22*:270, 1701.] This general law states the rate of heat loss of a body is directly proportional to the difference between the body temperature $T B D$ and a surrounding, ambient temperature $T B D$ .

Applied as a boundary condition, $T B D$ is the fluid temperature at the boundary, and $T B D$ a temperature some distance beyond the solid boundary. A heat transfer coefficient $T B D$ , with SI units $T B D$ , provides the constant of proportionality such that

$ T B D $
Substituting Eq. (4.28) and rearranging gives an equation for the *fixed heat transfer* *coefficient* condition:
$ T B D $
The equation has the form of a Robin condition, Eq. (4.10), so can be implemented as described in Sec. 4.9. The coefficient $T B D$ is typically characterised for the particular flow regime and solid boundary, by some estimate, experimental measurements or computer simulation.


== Recommended boundary conditions

This chapter covers a range of boundary conditions and their implementations. It first describes a specification of the basic conditions at inlet, outlet and wall boundaries for _subsonic_ flow with *fixed value* and *zero gradient*.

The conditions, based on the propagation of disturbances, are described in Sec. 4.3:

- zero gradient on $T B D$ at an inlet, fixed value on other variables;
- fixed value on $T B D$ at an outlet, zero gradient on other variables.


The conditions at a wall are similar to an inlet for $T B D$ and $T B D$ , but generally are represented more directly by physical models, _e.g._ the condition for heat flux $T B D$ for $T B D$ .

=== Supersonic conditions

The basic conditions for _supersonic_ flow are discussed in Sec. 4.5. If the flow speed is supersonic at an inlet, the basic condition is fixed value for $T B D$ ; it is zero gradient for $T B D$ if the flow is supersonic at an outlet.

=== Robust, practical conditions

Section 4.6 introduced a free boundary that cannot be defined as an inlet or outlet, but instead often uses the following conditions:

- *total pressure* for $T B D$ , see Sec. 4.7;
- *inlet-outlet-velocity* for $T B D$ , see Sec. 4.15;
- *inlet-outlet* for $T B D$ , see Sec. 4.10.


These conditions also respond well at an outlet, in the event that some inflow occurs at startup, a rotating structure passes through the boundary _etc._, see Sec. 4.10.

The *freestream* conditions, Sec. 4.16, are effective for cases with known $T B D$ and $T B D$ at a free, far-field boundary.

The *symmetry* and *wedge* conditions enable suitable cases to be simplified as symmetric and axisymmetric, respectively.

In the presence of a body force $T B D$ , the zero gradient condition for $T B D$ at inlets and walls should be replaced by a *fixed gradient* condition $T B D$ , see Sec. 4.4.

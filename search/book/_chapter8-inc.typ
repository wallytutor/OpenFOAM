The preceding chapters establish some general principles of CFD, covering topics in fluid dynamics, numerical methods, boundary conditions, algorithms, solvers and turbulence and its modelling.

This chapter aims to bring together some of these principles using four sample CFD problems, some of which are mentioned in the earlier chapters. The four examples demonstrate a variety of flow types, including:

- two _internal_ flows, inside a room and a Venturi tube, in which, beyond regions of significant inflow and outflow, the fluid is contained by solid boundaries;
- two _external_ flows, around a vehicle and a cylinder, in which a considerable part of the domain is unbounded.


The flow in the room and around the vehicle are turbulent, so demonstrate turbulence modelling and wall functions from Chapter 7. The flow is laminar in the other examples.

The vehicle and Venturi tube cases are simulated using the steady-state algorithm in Sec. 5.12, whereas the room and cylinder cases use the transient algorithm in in Sec. 5.19.



=== Mesh generation
The sample problems are supported by additional information on mesh generation. Strategies for meshing between the far-field and a body of interest in external flows are discussed in Sec. 8.1.
Meshes may include cells with high aspect ratio to increase the solution efficiency, particularly with external flows. Including these cells, however, increases the range of magnitudes of the matrix coefficients, which adversely affects solution convergence, as discussed in Sec. 8.2.

The generation of block-structured and unstructured meshes are discussed in Sec. 8.3 and Sec. 8.5 respectively. The block-structured approach is suited to problems whose boundaries are described by relatively simple geometric forms, _e.g._ cylinders, cuboids, _etc._

The unstructured approach is designed to accommodate more complex boundary geometries, _e.g._ the vehicle in the external aerodynamics example. The geometries are typically defined by a computer-aided design (CAD) model.



=== Boundary Conditions
The examples demonstrate a variety of specialised boundary conditions. The cylinder and vehicle both use the freestream condition, described in Sec. 4.16 at the far-field boundary. The Venturi example uses a nonuniform inlet condition for $T B D$ , described in Sec. 8.7. Finally, the flow in a room uses a condition for $T B D$ to ensure correct boundary fluxes, as described in Sec. 5.11 and the external heat flux condition for $T B D$ from Sec. 4.17.


== External flows

This chapter includes examples of flow around a cylinder and automotive aerodynamics. They are examples of _external_ flows, which are characterised by a lack of solid boundaries to define the extent of some part, or all, of the solution domain.

The cost of a CFD simulation increases with the size of the computational mesh, so we cannot afford to have a solution domain which extends indefinitely in some direction. Boundaries must to be drawn somewhere, but which locations are optimal?

To calculated forces on a body accurately, the external boundaries need to positioned a surprisingly long distance from the body. For example, CFD calculations of forces on an aerofoil have used a domain length of 1000 $T B D$ the aerofoil chord length.

To understand why large distances are necessary, we must first recognise that the pressure is imposed at external boundaries according to the applied boundary conditions. Thy can include a fixed value outflow condition in Sec. 4.3, a total pressure condition in Sec. 4.7, a freestream condition in Sec. 4.16, _etc._

The pressure at external boundaries will inevitably differ from the pressure at the same location in an unbounded problem. As discussed in Sec. 2.22, a change in pressure at any point influences the pressure everywhere in the solution domain, so any discrepancies in boundary pressures are reflected in the force calculations, including the shear forces.

Larger solution domains yield smaller discrepancies in pressure, producing solutions that better represent an unbounded domain. Accepting large domains are necessary, we need to limit the computational cost with an efficient meshing strategy that combines smaller cells to resolve higher gradients in complex, transient flow regions, with larger cells where the flow is more uniform.

=== Mesh refinement

An external flow approaches far-field freestream conditions in regions around the external boundaries. Larger cells are effective in these regions since the flow is fairly uniform; smaller cells can then be reserved for more critical regions, close to, and in the wake of, solid bodies.

Two common ways to vary cell size are shown above. To the left is an example of a structured mesh which increases the cell lengths towards the freestream boundaries. The transition in size is smooth but cells tend to be generated which are noticeably longer in one direction than another, _i.e._ with a high _aspect ratio_, discussed in more detail in the following section.

To the right is an example of a mesh which refines the mesh away from the freestream boundaries using cell splitting. Starting from a coarse mesh, the method reduces the cell size by splitting each cell into two in each direction. This approach produces more abrupt changes in cell size across refinement interfaces, but maintains lower aspect ratios in the cells.



== Aspect ratio

The previous section mentioned cells of a mesh that are longer in one direction than another. The parameter used to describe this phenomenon is _aspect ratio_, defined as the ratio of the longest to shortest dimensions, $T B D$ , of a cell.

Meshes include cells with high aspect ratio for efficiency, by exploiting a prior understanding that the gradients in some flow properties _e.g._ shear stress, deformation rate and velocity, are much greater in one direction than others.

They are most commonly used to resolve the boundary layers at walls. Greater mesh resolution is required normal to the wall, where gradients are particularly high, than tangential to the wall. The centre heights of cells adjacent to the boundary need to be very small, corresponding to $T B D$ .

In the aerofoil above, the aspect ratio is 1000 at the trailing edge. If those cells were replaced by ones with an aspect ratio of 1, with the same cell height, the number of cells would increase by a factor of 1000, illustrating the extent of the gains in efficiency.

Introducing high aspect ratio in one area of the mesh causes problems elsewhere, however. The aerofoil example maintains a high aspect ratio beyond the trailing edge, actually increasing to $T B D$ further downstream. If the flow deviates from a direction parallel to the long side of the high aspect ratio cells, transient simulations then require very small $T B D$ to fall within reasonable $T B D$ limits. The mesh can be modified by abruptly reducing the aspect ratio beyond the trailing edge, but this increases the error associated with discretisation.

With Laplacian discretisation by Eq. (3.2) and Eq. (3.5), the matrix coefficients include a factor $T B D$ . On a regular mesh, _both_ $T B D$ and $T B D$ are larger for the faces across the longest side of a high aspect ratio cell, than the shortest side, by a factor equal to the aspect ratio. Matrix coefficients can thereby differ by a factor of $T B D$ .

The convergence of the pressure equation, dominated by the Laplacian derivative, is adversely affected, in particular with the CG-based matrix solvers. The GAMG solver fares much better because the agglomeration strategy is based on a progressive reduction in aspect ratio as discussed in Sec. 5.18.

A matrix with a wide range of coefficient values yields lower residual values $T B D$ , calculated from Eq. (5.11). In other words, for a given $T B D$ , a solution with high aspect ratio cells is less converged than one without. Therefore, to maintain an equivalent level of convergence with high aspect ratio cells, tolerance controls need to be reduced, often by a few orders of magnitude.



== Block-structured meshes

This section describes the generation of meshes for cases with boundaries formed by relatively simple geometric shapes, _e.g._ a cylinder, such as the examples in Sec. 8.4 and Sec. 8.8.

Block-structured mesh generation involves splitting the domain (assumed 3D) into several 6-sided blocks. Each block is then divided into a number of hexahedral-shaped cells.

There are some popular strategies for specifying the structure of blocks around boundaries of particular shapes. The most general approach uses blocks which extend outwards from sections of the boundary into the domain. For a cylinder in cross-section, 4 blocks are typically used as shown below by the blocks whose edges are oriented at $T B D$ to the horizontal.

The domain extends away from the cylinder with additional blocks connected to those that encircle the cylinder. The mesh is usually graded starting with small cell heights on the cylinder surface, which gradually increase with distance from the surface. The main weakness of this block structure is the abrupt transition in non-orthogonality from $T B D$ to $T B D$ at the outer vertices of the inner blocks, indicated in the figure.

The figure above shows the structure of the mesh of 40,000 cells used in Sec. 8.4 for a cylinder of diameter $T B D$ . Each element represents 10 $T B D$ 10 = 100 cells in the actual mesh. There are 80 layers of cells around the cylinder, with grading that gives a centre height of $T B D$ in cells adjacent to the cylinder surface.

A mesh _inside_ a cylinder, _e.g._ a pipe, can use a similar structure in which the inner vertices of the curved blocks are connected to form a single block along the centre, as shown above.



== Flow around a cylinder

Flow around a cylinder was used as an example of boundary layer separation in Sec. 6.5. It showed an image of vortices shedding from the cylinder in a periodic manner, known as the Kármán vortex street.

The image comes from a CFD simulation in two dimensions, representative of a cylinder of infinite axial length, using the principal parameters: cylinder diameter $T B D$ ; freestream velocity $T B D$ in the $T B D$ -direction; and, fluid kinematic viscosity $T B D$ . The corresponding Reynolds number $T B D$ , which falls within the laminar flow regime.

The computational mesh used in the simulation is described in Sec. 8.3. The centre height of the cells adjacent to the cylinder was $T B D$ , corresponding to a calculated $T B D$ .

The simulation used the transient solution algorithm in Sec. 5.19, solving for momentum conservation for an incompressible fluid, with $T B D$ and constant $T B D$ . No energy equation was solved and no turbulence modelling was required since the flow was laminar.

The freestream boundary conditions from Sec. 4.16 were applied to $T B D$ and $T B D$ over the entire external boundary, with reference values $T B D$ , $T B D$ and $T B D$ . The no-slip condition, $T B D$ and $T B D$ , was applied on the cylinder boundary.

The simulation ran for $T B D$ with a time step $T B D$ using recommended numerical schemes from Sec. 3.23. Oscillations began in the wake of the cylinder at $T B D$ , soon leading to shedding of vortices which reached a stable pattern at $T B D$ .

The following figure, shaded by $T B D$ , shows the slower-moving vortices as the darker structures which propagate downstream of the cylinder.

The repeated vortex shedding causes oscillations in the force $T B D$ of the fluid _on_ the cylinder. The force is calculated as the sum of viscous and pressure forces acting on the cylinder patch faces ( $T B D$ ), in kinematic units $T B D$ (from kinematic $T B D$ and $T B D$ ) by

$ T B D $
The dimensionless lift coefficient $T B D$ is shown below, calculated from: the $T B D$ -component of the force is $T B D$ , where $T B D$ is the unit vector in the $T B D$ -direction; and, the projected frontal area $T B D$ of the cylinder.
The oscillation period $T B D$ corresponds to a Strouhal number $T B D$ , consistent with published data.#footnote[Christoffer Norberg, _Flow around a circular cylinder: Aspects of fluctuating lift_, 2001.]



== Unstructured hex-dominant meshes

While block-structured meshes are effective for model problems with simple boundaries, _e.g._ a cylinder, they are impractical for many engineering applications, _e.g._ automotive aerodynamics.

Instead, CFD practitioners turn to tools which can generate meshes of many millions of cells and accommodate complex boundaries, _e.g._ a from a CAD model of a vehicle, in a highly automated manner.

A logical approach is to apply hexahedral-shaped cells that are oriented with the known, or anticipated, direction of flow. This is possible towards far-field boundaries and along solid walls, with the far-field cells generally larger than those at solid walls.

Cells must then be generated to fill the region between the aligned meshes near the far-field and solid boundaries. One approach is to infill the region with tetrahedral-shaped cells using automated methods based on Delaunay triangulation.#footnote[Boris Delaunay, _Sur la sph__ére vide_, 1934.]

A more recent strategy is to fill the entire domain with large cells, similar to those in the far-field, which are then progressively split until they reach the desired size at solid boundaries. Subsequent steps then conform and align the cells to the solid boundaries.

The infill region then contains mainly hexahedra, with the mesh described as “hex-dominant”. Where there is transition in levels of cell splitting (see above), the larger cells become polyhedrons with an increased number of faces due to the smaller cells.

The splitting reduces solution accuracy by introducing a modest amount of non-orthogonality at the intersecting faces, _i.e._ $T B D$ for cells of unit aspect ratio. It also causes the weights $T B D$ for interpolation schemes, see Sec. 3.7, to deviate from the ideal $T B D$ ; instead, for cells of unit aspect ratio, $T B D$ (or $T B D$ ).

The generation of unstructured hex-dominant meshes is highly automatable, using simple controls that define the levels of cell splitting within specified regions of the domain, as shown above.



== Automotive aerodynamics

An example of flow around a road vehicle was used to discuss some boundary conditions in Chapter 4 and to illustrate the cost of simulating turbulence in Sec. 6.8. An _aerodynamics_ simulation was undertaken to capture the air flow around the vehicle, described by a CAD model. The aim was to calculate the drag coefficient at a speed of $T B D$ .

A mesh of 20 million cells was generated, with the vehicle facing a freestream flow velocity $T B D$ . The vehicle and ground formed solid boundaries, with far-field boundaries positioned $T B D$ upstream and $T B D$ downstream of the vehicle.

Along the elevated sections of the far-field boundary, the cell length was $T B D$ , reducing to $T B D$ towards the vehicle by splitting within specified regions. Additional cell layers along the vehicle surface resulted in a near-wall cell height of $T B D$ .

The simulation used the steady-state algorithm in Sec. 5.12, with an incompressible fluid with uniform $T B D$ .

The freestream boundary conditions from Sec. 4.16 were applied to $T B D$ and $T B D$ at the far-field boundaries, with reference values $T B D$ , $T B D$ and $T B D$ . The condition $T B D$ was a applied at solid boundaries, with $T B D$ applied to the vehicle and $T B D$ on the ground to emulate their relative motion.

Turbulence was modelled using the $T B D$ SST model described in Sec. 7.11. Turbulence levels of $T B D$ and $T B D$ were applied at the freestream boundaries and the standard wall function from Sec. 7.5 was applied at the vehicle and ground.

The simulation ran for 3000 iterations using numerical schemes recommended in Sec. 3.23. The drag coefficient $T B D$ was calculated from the projected frontal area $T B D$ and the $T B D$ -component $T B D$ of the force $T B D$ on the vehicle using Eq. (8.1).

The flow in the wake of the vehicle is naturally unsteady, which prevents convergence to a steady-state solution. Beyond 1500 iterations, however, the solution oscillates around an estimated mean $T B D$ .



== Nonuniform inlet velocity

At an inlet boundary, the velocity is usually specified by a fixed value condition $T B D$ , as discussed in Sec. 4.3. In many CFD applications, the flow at an inlet boundary is described by a single speed which is applied to all faces assuming $T B D$ is uniform across the boundary.

This creates an anomaly where the inlet boundary meets a wall. In the vicinity of the two boundaries, the flow decelerates from $T B D$ at an inlet face to a value in the adjacent cell close to the no-slip condition $T B D$ applied at the wall.

There is inevitably a high “spike” in pressure and shear stress within the cell in order to decelerate the flow so rapidly. As the length of the cell (in the flow direction) is reduced, the deceleration and associated pressure $T B D$ increases such that $T B D$ in the limit that cell volume $T B D$ .

The solution tends to converge more slowly with the pressure spike, and can be unstable. Furthermore, the spike in shear stress can generate high levels of turbulence which can cause flow separation where the inlet and wall meet.

The uniform condition does not reflect the flow behaviour upstream of the inlet. For example, assuming the wall extends upstream, a boundary layer would have developed at the inlet; or if the wall begins at the inlet, the flow would stagnate at its leading edge.

This is a good example of the axiom that numerical methods do not respond well to any modelling which is unphysical. The problems can be avoided by specifying a nonuniform $T B D$ which represents the upstream flow better.

Some fields of engineering use established theories to describe the inlet $T B D$ , _e.g._ wind engineering uses a profile of $T B D$ based on an atmospheric boundary layer along the earth’s surface.

More generally, a nonuniform $T B D$ can be specified which tends to $T B D$ at wall boundaries. Profiles can be described by $T B D$ where: $T B D$ is the normalised velocity and $T B D$ is normalised distance to the nearest wall boundary; and, $T B D$ and $T B D$ denote the maximum $T B D$ and $T B D$ values.

It is logical to specify $T B D$ using established profiles for boundary layers as a reasonable estimate for the upstream conditions. A quadratic function $T B D$ represents a developed boundary layer for laminar flow, matching the analytical profiles for flow in a pipe or between flat parallel plates, _i.e._ Poiseuille’s law.

Alternatively, a power law function $T B D$ represents a developed turbulent boundary layer quite well. Prandtl used $T B D$ — his _one-seventh power_ _law_#footnote[Ludwig Prandtl, _The mechanics of viscous fluids_, 1935.] — to reproduce data for flow in a pipe, but any suitable exponent $T B D$ can be used in practice.



== Venturi tube

A Venturi tube is a device used to measure the flow rate of a fluid in a pipe. The device exploits the Venturi effect#footnote[Giovanni Battista Venturi, _Recherches exp__érimentales sur le principe de la communication_ _lat__érale du mouvement dans les fluides appliqu__é a l’explication de diff__érens ph__énom__ènes hydrauliques_, 1797.], _i.e._ that pressure reduces when the fluid flows through a restriction.

The volumetric flow rate $T B D$ is calculated according to

$ T B D $
where: $T B D$ is the decrease in _kinematic_ pressure $T B D$ between the Venturi inlet and throat; $T B D$ and $T B D$ are the inlet and throat diameter, respectively; the cross-sectional area is $T B D$ ; and, $T B D$ is a _discharge coefficient_. A coefficient of $T B D$ satisfies the Bernoulli equation, $T B D$ , between the inlet and throat, where $T B D$ is the fluid speed.
For $T B D$ , $T B D$ is calculated accurately using $T B D$ values between 0.97 and 1,#footnote[ISO 5167-4, _Measurement of fluid flow by means of pressure differential devices inserted in circular_ _cross-section conduits running full – part 4: venturi tubes_, 2003.] but, for $T B D$ , suitable values of $T B D$ decrease significantly below 1 in order to account for pressure losses due to viscous forces.

A simulation was performed to calculate $T B D$ for a Venturi tube, shown in the figure, with $T B D$ and an inlet $T B D$ . The inlet velocity was specified using the quadratic profile in Sec. 8.7 with a mean cross-sectional speed $T B D$ , corresponding to a maximum speed $T B D$ .

The simulation used the steady-state algorithm in Sec. 5.12, with an incompressible fluid with uniform $T B D$ . The mesh contained 57,600 cells, with a near-wall cell height of 1.5mm, which resolved the velocity profiles as shown below. The flow recirculates near the wall downstream of the Venturi throat, causing inflow at the outlet boundary. Consequently, the total pressure and inlet-outlet-velocity boundary conditions, described in Sec. 4.7 and Sec. 4.15 respectively, were applied at the outlet to maintain stability.

The flow was laminar so no turbulence modelling was required. The solution converged to within an absolute tolerance $T B D$ , see Sec. 5.4, in 292 iterations.

The pressure drop was $T B D$ between the centres of the inlet and throat sections, with a corresponding $T B D$ .



== Heating a room

A simulation of heating a room demonstrates natural convection driven by buoyancy forces. An idealised, ground floor room is presented below, with external glass doors and sloping roof, and internal walls and ceiling.

A heater is located along one side wall, below the point where the roof and ceiling meet. The aim of the simulation was to calculate the room temperature with the heater running at 1kW, when the ambient external temperature $T B D$ .

The thermal boundary conditions were specified as follows: floor temperature $T B D$ ; ceiling $T B D$ , representing the first floor temperature; insulated walls, with $T B D$ , see Sec. 4.17; glass doors and roof with a heat flux according to Eq. (4.30) using $T B D$ and $T B D$ , respectively.

The mesh contained 350,000 hexahedral cells with grading that gave a cell height of approximately $T B D$ along the walls.

Transport properties for air, $T B D$ and $T B D$ , were used. Turbulence was modelled using the $T B D$ SST model described in Sec. 7.11, with initial levels of $T B D$ and $T B D$ . The near-wall cell centres corresponded to $T B D$ , so the continuous wall function from Sec. 7.6 and thermal wall function from Sec. 7.14 were applied at the boundaries.

The simulation used the transient solution algorithm in Sec. 5.19, including the buoyancy force $T B D$ in Eq. (2.67), with $T B D$ . The condition $T B D$ was applied at all boundaries, combined with the flux calculation in Eq. (5.20).

The variations in $T B D$ within $T B D$ were calculated using the ideal gas Eq. (2.55) using $T B D$ .

The simulation ran with a time step $T B D$ . The flow is highly unsteady, but at $T B D$ the heat losses through the boundaries oscillate about the mean levels indicated above.

Between the ground and 2m, occupants experience a variation in $T B D$ . In the space adjacent to the roof and ceiling, the higher $T B D$ generates significant heat losses, especially to the outside through the roof.



== Building a CFD simulation

The examples in this chapter aim to demonstrate the general principles of CFD described in the previous chapters. They also highlight other important aspects of CFD, _e.g._ mesh generation. To finish with, here is some advice on building CFD simulations.

Before doing a simulation it is important to define its purpose, including the information you want to obtain. Initial data about the problem should be gathered, including the solution domain and boundaries, flow conditions, and models and properties that might be used. The flow $T B D$ should be estimated to establish whether the flow is laminar or turbulent (Chapter 6).

It is important to recognise that effective CFD simulations are not created in one attempt. Instead, they follow a typical design process, as shown below.

The process is iterative, beginning with a _prototype_ simulation. The prototype should be as _simple_ as possible in order to reach the first successful test quickly. Once achieved, the design cycle is set in motion and, thereafter, the simulation can be evaluated and improved in incremental steps, following the cycle above. Any problem can be attributed to the most recent change, which is much easier to diagnose when the change is small.

Simulations need to run quickly so that frequent, small changes can be tested efficiently. They run in a few minutes on a mesh of $T B D$ cells, which is a good initial size for the prototype. Steady-state solutions (Chapter 5) run particularly quickly.

Fields, _e.g._ $T B D$ , must be initialised and the boundary conditions applied (Chapter 4) early in the design process. The prototype can include simple conditions, _e.g._ fixed value for $T B D$ , before switching to more complex conditions, _e.g._ a heat flux for $T B D$ .

If the flow is turbulent, robust RAS models should be deployed initially (Chapter 7). The prototype mesh size dictates that boundary layers are inevitably modelled using wall functions.

If a prototype simulation fails to run, numerical causes should be investigated, starting with mesh quality, _e.g._ non-orthogonality. Problems with discretisation schemes can be eliminated by applying the most stable schemes first, _e.g._ upwind (Chapter 3).

The initial simulations may simply establish a basic flow, solving mass and momentum conservation (Chapter 2) with simplifying assumptions, _e.g._ the incompressibility condition. The physical models should be simple, using constant properties.

The incremental changes to the simulation then incorporate additional models and equations. Often they require additional fields, boundary conditions, discretisation schemes, _etc._

Beyond that, the simulation may venture into more complex areas of CFD modelling, including multiphase flows, conjugate heat transfer, compressible flows, reactions, particle methods and large-eddy simulation — outside the scope of this book.

However, the general principles presented in this book still apply whether a CFD simulation is simple or complex. And when simulations are more complex, they fail just as often because of a lack of adherence to these principles as for any other reason.

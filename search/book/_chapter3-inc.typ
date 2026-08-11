Chapter 2 described:

- the main governing equations of fluid dynamics, namely conservation laws of mass, momentum and energy;
- standard models for diffusion of momentum and heat, namely Newton’s law of viscosity and Fourier’s law of conduction;
- the ideal gas equation of state and specific heat capacity;
- tensor algebra used in the equations and the derivatives of time, divergence, gradient and Laplacian.


Sets of governing equations and models are combined to simulate different flow problems. Numerical methods are used in CFD to solve the equations since closed-form, analytical solutions only exist for very simple geometries and flow regimes.

The numerical methods represent _continuous_ physical entities by equivalent _discrete_ entities, _i.e._:

- _time_ is split into intervals of duration $T B D$ ;
- _space_ becomes cells of a mesh — the solution domain;
- _fields_, _e.g._ $T B D$ , become discrete values, _e.g._ one per cell;
- _partial differential equations_, _e.g._ for momentum, become a set of linear equations.


There are numerous discretisation methods to create the set of linear equations, but this book specifically examines the finite volume method.

The sets of equations can be solved using a variety of algorithms, with finite volume discretisation lending itself towards iterative methods.



== The finite volume concept

There are many numerical methods which can be used to solve the partial differential equations encountered in fluid dynamics. No single method is unequivocally better than others. Rather, the efficacy of each method depends on aspects of the case being simulated, _e.g._ its size, the required level of accuracy and the characteristics and complexity of the equations being solved.

For a method to be useful, it needs to be programmed into software. Complex programs require good design and design relies on good concepts. So inevitably the concepts behind any method are as important as the details contained within.

The finite volume method adopts the idea of control volumes used to create models of physical systems. A control volume represents a region of space, which is generally fixed, enclosed by a surface through which fluid flows in and out.

It applies conservation equations, _e.g._ of mass, momentum and energy, by balancing fluxes, due to inflow and outflow at the bounding surface, with additional sources within the volume. Rather than using a single control volume to describe an entire physical system, _e.g._ a heating tank, the finite volume method splits the system domain, _i.e._ the tank, into multiple connected finite volumes. Conservation equations are applied to each volume, ensuring that the fluxes of mass, momentum and heat across surfaces are consistent between the volumes they connect.

The perceived wisdom is that the finite volume method was first introduced in the early 1970’s. But the PhD thesis of Runchal#footnote[Akshai Runchal, _Transfer processes in steady two-dimensional separated flows_, PhD Thesis, 1969.] from 1969 describes a method which is clearly finite volume. He also published the figure reproduced below,#footnote[Akshai Runchal, Brian Spalding and Micha Wolfshtein, _Numerical solution of the elliptic equations for_ _transport of vorticity, heat, and matter in two-dimensional flow_, 1969.] which displays the computational mesh as a set of connected control volumes.

Runchal credits the idea to his PhD supervisor Brian Spalding who provided an analogy of tanks connected by tubes in 1967. It captures the essence of finite volume, which uses a mesh to define a physical system of control volumes. Its emphasis is on calculating fluxes between volumes and ensuring conservation.

By contrast, other methods, _e.g._ finite element, use the mesh to construct mathematical functions to calculate distributions of properties. The finite volume method does not do that.



== Computational mesh

A CFD simulation begins with a _solution domain_ which specifies a region of space of a particular geometric shape, in which fluid dynamics equations are solved. For example, when simulating flow along a pipe, the solution domain would be some region of space occupied by the fluid within the pipe.

The process of _mesh generation_ subdivides the solution domain into a mesh of small volumes, or _cells_. Computer programs are used to create the cells according to some user specification. We discuss neither mesh generation software nor any underlying methods they use. Rather, we define the structure and properties of a valid mesh in this chapter and some aspects of the design of optimal meshes in Chapter 8.

=== Polyhedral mesh

Meshes in modern CFD software, that use the finite volume method, can contain cells of any irregular polyhedral shape. A cell can consist of an unlimited number of faces ( $T B D$  4) and each face can have an unlimited number of edges ( $T B D$  3).

Cells are _contiguous_, _i.e._ the faces of a given cell are common to its neighbouring cells (unless they form the boundary of the solution domain). There is no restriction on the alignment of cells with co-ordinate axes.

The figure above shows two cells connected by a common face displayed in grey. Every face is defined by a list of vertices that follow edges in sequence around the face.

Each face inside the mesh is common to two cells, one known as the _owner_ cell, the other the _neighbour_ cell. At the domain boundary, each face is connected to one owner cell. By assigning a unique index to each cell, the connectivity between cells can represented by storing the indices of the owner and neighbour cell of each face.

The cell faces that describe the domain boundary are separated into groups, each with a unique name. Each named group of boundary faces, known as a _patch_, identifies a specific region of the domain boundary, to enable us to apply specific boundary conditions to that region when running a CFD simulation.

A mesh is therefore defined by the following:

- a list of every vertex used to define all faces;
- a list of faces, each defined by a sequence of vertex indices;
- a list of owner cell indices associated with the faces;
- a list of neighbour cell indices associated with faces;
- grouping boundary faces into a set of patches.




== Finite volume mesh

The finite volume numerical method is closely associated with the concept of surface and volume integrals used in Chapter 2. Below, we extract the main geometric elements from the figures used for the conservation laws.

The integrals use volumes $T B D$ and area vectors $T B D$ . The numerical method uses equivalent discrete quantities for cells and faces:

- cell volume $T B D$ ;
- face area vector $T B D$ , with area magnitude $T B D$ ;
- face unit normal vector $T B D$ .


The finite volume method relates discrete values of fields, _e.g._ pressure, to cells and faces within the mesh. For many calculations data must to be assigned to point locations, in particular the cell centre (more specifically _centroid_) $T B D$ and face centre $T B D$ .

=== Calculating mesh data

To calculate $T B D$ , each polygonal face is decomposed into triangles using an apex point $T B D$ . The area vector $T B D$ and centre $T B D$ are then calculated for each triangle according to

$ T B D $
where “ $T B D$ ” is the cross product, Eq. (2.70). The sum of area vectors gives $T B D$ and $T B D$ is the area weighted sum of triangle centres over $T B D$ triangles, calculated by
$ T B D $
The cell volume $T B D$ can be calculated from $T B D$ and $T B D$ by Gauss’s theorem, using $T B D$ to describe position and noting $T B D$ . The surface integral ( $T B D$ ) becomes a sum over cell faces $T B D$ , replacing discrete values $T B D$ and $T B D$ for $T B D$ and $T B D$ respectively as follows:
$ T B D $
The cell centre is calculated similarly, noting $T B D$ , by
$ T B D $



== Equation discretisation

Equation discretisation converts partial differential equations for _continuous_ fields, _e.g._ pressure $T B D$ , into sets of linear equations for _discrete_ fields.

The values of the principal fields, _e.g._ $T B D$ , are associated with cells. A field is then represented by an array of values, $T B D$ , for cell indices $T B D$ . $T B D$ is the total number of cells.

Equation discretisation creates a linear equation for each cell. For cell 43 above, the equation might have the following form:

$ T B D $

where $T B D$ and $T B D$ are coefficients corresponding to cell indices $T B D$ , $T B D$ (diagonal coefficient in *bold*). The set of linear equations for all cells can be written as a matrix equation of the form:

$ T B D $
The matrix contains of a set of coefficients $T B D$ where each row $T B D$ corresponds to the linear equation for the cell with index $T B D$ . Each row of coefficients are non-zero only for the respective cell (diagonal $T B D$ ) and near-neighbours. All other coefficients are zero, making the matrix extremely _sparse_.
The matrix equation can be represented as $T B D$ , where: $T B D$ are the matrix coefficients; $T B D$ are the source coefficients $T B D$ ; and, $T B D$ is the discretised pressure field. It may also be illustrated as follows, where ‘ $T B D$ ’ indicates non-zero coefficients.

$ T B D $

=== Segregated solution

A CFD simulation generally solves a set of physical equations, _e.g._ for mass, momentum and energy conservation. The finite volume method traditionally discretises each physical equation separately to form individual matrix equations for single solution variables, _e.g._ $T B D$ , $T B D$ , $T B D$ , rather than creating a single matrix equation that represents all the physical equations.

The _segregated_ matrix equations are solved one variable at a time, _e.g._ solving for $T B D$ , $T B D$ and $T B D$ in separate steps. Where the solution variable is a vector or tensor, _e.g._ $T B D$ , it is _decoupled_ into individual matrix equations for each component, _e.g._ $T B D$ , $T B D$ , $T B D$ .

The matrix equations are solved in an _iterative_ sequence, in which the equation for one variable, _e.g._ $T B D$ , incorporates current values of other variables, _e.g._ $T B D$ , $T B D$ and $T B D$ , into the source vector, as shown below.



== Matrix construction

The construction of each matrix equation involves building the matrix and source coefficients from the terms in the equation being solved, with further adjustments for boundary conditions.

The figure below illustrates the process of building a matrix equation for a field $T B D$ from an equation including advection, diffusion and source $T B D$ .


Coefficients for matrix $T B D$ and source $T B D$ are calculated for each individual term in the equation, _e.g._ $T B D$ , $T B D$ , _etc._ using the discretisation methods described in this chapter.
The overall coefficients are calculated as the sum of coefficients for each term in the the equation. Most terms contribute to both the matrix and source coefficients, although this depends on the choice of discretisation scheme.

Finally, boundary conditions are incorporated into the equation through further adjustments to coefficients in $T B D$ and $T B D$ as shown below. The adjustments, principally from the advection and diffusion terms, are applied to coefficients corresponding to cells at the domain boundary.



=== Implicit and explicit

The equation for a field $T B D$ on the previous page is discretised to form the matrix equation $T B D$ . The discretisation of a term is _implicit_ when it contributes to coefficients in $T B D$ by treating $T B D$ as the solved field $T B D$ .

_Explicit_ discretisation calculates coefficients in $T B D$ only, using current values of fields. When solving an equation for $T B D$ , derivatives without $T B D$ must be explicit. Terms with $T B D$ could be treated explicitly by using current values of $T B D$ , but they generally are not, since explicit solutions are unstable beyond a limiting time step as described in Sec. 3.17. A notable exception to this are the terms discussed in Sec. 3.20.

The curl derivative, _e.g._ $T B D$ , includes terms in $T B D$ and $T B D$ in the decoupled matrix equation for $T B D$ . These terms must be treated explicitly since they do not include $T B D$ itself. The situation for $T B D$ and $T B D$ is the same, so the curl derivative can be explicit only.

This leaves the following terms which are generally treated implicitly:#footnote[also equivalent terms including density, _e.g._ $T B D$ , $T B D$]

- time derivative $T B D$ ;
- diffusion (Laplacian) $T B D$ ;
- advection $T B D$ ;
- implicit linear function $T B D$ , where $T B D$ is a scalar.


This chapter details: the discretisation of these terms, which are generally treated implicitly; and, other terms such as divergence $T B D$ and gradient $T B D$ which can only be discretised explicitly within a segregated solution.



== Overview of discretisation

The choice of numerical method determines how coefficients for $T B D$ and $T B D$ are calculated and consequently the characteristics of the resulting matrix equation $T B D$ .

The finite volume method described here is firmly rooted in the underlying concept of control volumes, described in Sec. 3.1. It uses integrals over a surface surrounding a volume applied to irregular polyhedral meshes, described in Sec. 3.2.

The discretisation is described in terms of differential operators, _e.g._ $T B D$ and $T B D$ , applied to a general field $T B D$ .#footnote[consistent with the idea of “Field Operation and Manipulation”, which forms the acronym “FOAM” (later, OpenFOAM) used for the CFD software created by author Henry Weller.]

The main concept is that faces within a mesh form closed surfaces surrounding finite volumes, _e.g._ a single cell. Any surface integral that represents a derivative, _e.g._ $T B D$ , is approximated by a summation over the faces ‘ $T B D$ ’ that form the surface, _i.e._

$ T B D $
The flux associated with $T B D$ , _i.e._ $T B D$ , must then be calculated. The value $T B D$ is required at each _face_, which must be calculated by some method of interpolation of values of $T B D$ from cells neighbouring the respective face.
=== Intensive and extensive properties

In this chapter, derivatives and their discretisation are described _at a point_, _e.g._ $T B D$ , _e.g._ Eq. (3.8):

$ T B D $
Calculating a derivative with this expression using the known field $T B D$ simply produces another field (with values defined at cell centres, with units of $T B D$ /time).
The resulting field is $T B D$ , meaning it is independent of the size of the system/geometry. Like other intensive fields, _e.g._ $T B D$ and $T B D$ themselves, it can be used in further calculations, _e.g._ within another derivative or added/subtracted from other fields.

_Extensive_ properties are dependent on the size of the system. For example, the volumetric flux $T B D$ described in Sec. 3.9 is dependent on the face areas $T B D$ . Numerical operations involving extensive properties, _e.g._ addition, subtraction, or mapping to another location, generally produce meaningless data.

While _calculations_ of derivatives yield intensive fields, _matrix equations are constructed in_ _extensive form_, with coefficients and source vector scaled by cell volume $T B D$ . In other words, in the discretisation example above, the multiplication by $T B D$ would be omitted.

There is no $T B D$ multiplier in the discretisation of terms which do not involve a surface integral, _e.g._ time derivative Eq. (3.21) and terms in Sec. 3.20. For those terms, the calculation of matrix coefficients and sources includes a multiplication by $T B D$ .



== Laplacian discretisation

Let us first describe the discretisation of the Laplacian term for diffusion, introduced in Sec. 2.14.

Following the finite volume principles described in Sec. 3.1, the discretisation approximates the surface integral by a summation over faces by

$ T B D $

The $T B D$ in Eq. (3.2) can be viewed as the summation over faces of a single cell. Applying the summation to all cells provides the contribution to coefficients $T B D$ and $T B D$ of a matrix equation.

The mesh data $T B D$ and $T B D$ are calculated according to Sec. 3.3 so the remaining properties to be determined are:

- the diffusivity at faces $T B D$ ;
- the surface normal gradient at faces $T B D$ .


Fields $T B D$ and $T B D$ are associated with cells, so numerical schemes are required to evaluate properties at faces. We will first describe interpolation for $T B D$ , for which the linear scheme is generally used. The surface normal gradient $T B D$ is discussed in Sec. 3.8.

=== Interpolation from cells to faces

Mapping data between different locations is a common practice in numerics. Since the finite volume method is concerned with fluxes at faces, the principal mapping procedure is _interpolation_ from cells to faces. Since other interpolations are much less common, it can be assumed that the term “interpolation” means from cell to face unless stated otherwise.

For irregular polyhedral meshes, interpolation is generalised by defining a _weights_ field $T B D$ for each face according to

$ T B D $
where $T B D$ is the interpolated face field. The subscripts $T B D$ and $T B D$ indicate values at owner and neighbour cells, respectively.
=== Linear interpolation

The *linear* interpolation scheme sets $T B D$ according to a linear variation between cells values $T B D$ and $T B D$ . The weights can then be calculated based on distances from the face centre to adjacent cell centres, in the direction normal to the face, by

$ T B D $



== Surface normal gradient

The surface normal gradient $T B D$ is a part of the Laplacian discretisation Eq. (3.2), illustrated in the figure below.

The discretisation of $T B D$ is built upon a finite difference of cell values on each side of the face according to

$ T B D $
where $T B D$ . When this *orthogonal* scheme is applied to Eq. (3.2) to discretise a Laplacian, it forms coefficients $T B D$ of a matrix equation $T B D$ since it references cell values of the field $T B D$ . For cell $T B D$ , the coefficient for each neighbour cell ( $T B D$ ) is $T B D$ and the diagonal coefficient is the negative of the sum of neighbour coefficients: $T B D$ .
Discretisation of $T B D$ by Eq. (3.5) is most accurate when the face is orthogonal to $T B D$ , _i.e._ the angle $T B D$ between $T B D$ and $T B D$ is zero. However, if the face is _non-orthogonal_, the error associated with Eq. (3.5) increases with $T B D$ .

=== Non-orthogonal correction

A more accurate discretisation of $T B D$ at a non-orthogonal face is formed of the vector sum of the orthogonal scheme $T B D$ and an _explicit_ correction $T B D$ . The latter is calculated from the full gradient $T B D$ in adjacent cells (described in Sec. 3.15), interpolated to the face $T B D$ .

The correction $T B D$ is explicit, _i.e._ calculated using known values of $T B D$ , so may need updating within an iterative sequence to maintain accuracy, as discussed in Sec. 5.20. To ensure that the iterative sequence converges, the implicit contribution is elevated by replacing $T B D$ in the orthogonal scheme with

$ T B D $
The *corrected* $T B D$ scheme combines the implicit and explicit parts by
$ T B D $
The corrected scheme is generally stable for $T B D$ . For $T B D$ , stability can be maintained at the expense of accuracy by _limiting_ the magnitude of the correction $T B D$ below some fraction of the magnitude of the orthogonal $T B D$ part.


== Advection discretisation

In Sec. 2.8, we described advection terms as those of the form $T B D$ or $T B D$ . The inclusion of velocity $T B D$ within the divergence gives the term particular characteristics that require special treatment in discretisation.

Following the finite volume principles outlined in Sec. 3.1, the discretisation approximates the surface integral by a summation over faces by

$ T B D $

The discretisation requires calculation of the volumetric flux $T B D$ (see Sec. 2.3). In the case of $T B D$ included in the advection term $T B D$ , $T B D$ is the mass flux $T B D$ .

In the flux calculation, the interpolation of $T B D$ at cell centres to $T B D$ at faces uses the linear scheme of Sec. 3.7. Similarly, the linear scheme is used for the interpolation $T B D$ .

The critical issue — one of the most important in CFD numerics — is how to express our advected property $T B D$ at a face in terms of values $T B D$ in neighbouring cells.

=== Advection scheme introduction

The advected property $T B D$ is transported in the direction of flow velocity $T B D$ . Interpolation to the face of $T B D$ usually involves the flow direction. In the graphic below a face f is positioned between two cells. Based on the flow direction, the cells are labelled upwind ‘U’ and downwind ‘D’.

The *linear* interpolation scheme, #footnote[also known as the _central difference_ scheme, particularly in the context of advection.] described in Sec. 3.7, does not use the flow direction but expresses $T B D$ in terms of adjacent cells. At first sight, this choice of scheme is logical when considering _accuracy_. However, for advection, the linear scheme tends to generate unbounded solutions which are unstable.

The *upwind* scheme simply represents $T B D$ by the value in the upwind cell $T B D$ . It makes sense from a physical perspective since particles of fluid in the upwind cell are destined to travel to the face, transporting property $T B D$ with them.#footnote[or, to quote Brian Spalding: “the wind from a pigsty always stinks”.]

While the linear scheme is generally unbounded, the upwind scheme exhibits poor accuracy. In following section we explore the behaviour of upwind before looking at schemes that offer greater accuracy while attempting to maintain boundedness.



== Upwind scheme

The upwind scheme represents the face value $T B D$ by the value $T B D$ in the cell upwind of the face. The advantage of upwind is that it can _guarantee_ boundedness of a field $T B D$ . We can demonstrate this point by revisiting the 1D Eq. (2.32) in Sec. 2.9:

$ T B D $
In the graphic above, we track the translation of a profile of $T B D$ by equating changes in $T B D$ in time to the local gradient $T B D$ . If we apply upwind to calculate the change at point P, the gradient $T B D$ and no change in $T B D$ is correctly calculated (left).
However, linear differencing between upwind and downwind values results in $T B D$ , so predicts a decrease in the value at P (right). The solution produces a solution with $T B D$ , so is unbounded.

Boundedness of the conservative form of advection $T B D$ is only guaranteed when $T B D$ , as discussed in Sec. 2.9. In 1D, the conservative form moves $T B D$ inside the derivative $T B D$ . That gradient is only zero with upwind when $T B D$ is uniform, _i.e._ the 1D equivalent to $T B D$ .

=== Diffusion of upwind

The upwind scheme is highly _diffusive_ which can result in poor accuracy. Its diffusive nature can be explained by considering the following Taylor’s series expansion:

$ T B D $
In our 1D example, the upwind scheme calculates $T B D$ using $T B D$ and $T B D$ at locations U and P, separated by distance $T B D$ . Relating the upwind calculation to Eq. (3.9) gives
$ T B D $
In other words, the upwind discretisation represents $T B D$ but also the second derivative $T B D$ (and higher derivatives). $T B D$ is equivalent to a Laplacian, described in Sec. 2.14, which diffuses $T B D$ with a diffusivity proportional to $T B D$ .
The upwind scheme is particularly diffusive when the flow direction is not aligned with the cells of a mesh. In the 2D box of cells above, $T B D$ is advected at a $T B D$ angle, beginning with an abrupt step change from $T B D$  = 1 and $T B D$  = 0 between the left and lower boundaries. The step rapidly diffuses along the direction of travel as shown in graph (right) and shaded area (left).



== Limited advection schemes

Alternative schemes for advection attempt to overcome problems with boundedness and accuracy of the linear and upwind schemes respectively. Many schemes apply a _limiter_ $T B D$ between $T B D$ from upwind and $T B D$ from the linear scheme Eq. (3.4) according to

$ T B D $
When $T B D$ , the scheme reduces to upwind and it becomes linear interpolation when $T B D$ . For a uniform mesh ( $T B D$ ), $T B D$ represents interpolation using the downwind cell value.
Limited schemes attempt to optimise $T B D$ at each face, based on the local $T B D$ , to maximise accuracy whilst maintaining boundedness.

Many schemes analyse the change in gradient of $T B D$ between the face and upwind cell in the direction $T B D$ connecting cell centres. They define a function of a ratio $T B D$ of consecutive gradients as:

$ T B D $
There are numerous published schemes that define the limiter $T B D$ as a function of the gradient ratio $T B D$ . Those schemes that are most useful are described in the following sections.
=== Total variation diminishing schemes

Many useful schemes fall into a class known as Total Variation Diminishing (TVD).#footnote[Ami Harten, _High resolution schemes for hyperbolic conservation laws_, 1983.] The TVD idea is that if the total variation of field $T B D$ does not increase in time, “overshoots” and oscillations associated with unboundedness will not occur.

To qualify as TVD, the limiter function $T B D$ must fall within the shaded area in a Sweby diagram (above).#footnote[Peter Sweby, _High resolution schemes using flux limiters for hyperbolic conservation laws_, 1984.] The TVD concept is a 1D analysis. For 3D CFD on irregular polyhedral meshes, oscillations are more likely to occur with TVD schemes whose $T B D$ functions tend significantly to downwind, _i.e._ towards the upper part of the shaded area near $T B D$ .

A further property of a limited scheme is _symmetry_. A scheme is symmetric when the condition $T B D$ is satisfied. When this occurs, the scheme applies the same limiter to the gradient of $T B D$ , irrespective of the _sign_ of its gradient. As a consequence, a property $T B D$ , initialised with a symmetric profile, _e.g._ a bell curve, will retain its symmetry under advection.



== Useful TVD schemes

For irregular 3D polyhedral meshes, only TVD schemes whose limiter function $T B D$ avoids excessive downwind interpolation, are sufficiently robust to be practically useful.

=== Limited linear and minmod

The *limited linear* scheme — as the name suggests — maintains the linear scheme as much as possible to be TVD, before limiting to upwind to improve boundedness, according to

$ T B D $
Limited linear is a robust scheme that exhibits low diffusion, but has a disadvantage that it is not symmetric, according to the criterion described in Sec. 3.11.
Like limited linear, the *minmod* scheme#footnote[Philip Roe, _Characteristic-based schemes for the Euler equations_, 1986.] does not introduce any downwind interpolation. It is symmetric but is more diffusive by limiting more readily towards upwind, according to

$ T B D $

=== van Leer and van Albada

There are two notable TVD schemes that: exhibit low diffusion; maintain symmetry through some downwind interpolation; and, have limiter functions which are continuous in gradient which avoids problems due to sudden changes in limiter values.

The *van Leer* scheme#footnote[Bram van Leer, _Towards the ultimate conservative difference scheme II. Monotonicity and conservation_ _combined in a second order scheme_, 1974.] is the original limited advection scheme described by the limiter function

$ T B D $
Since $T B D$ as $T B D$ , the van Leer scheme uses relatively strong downwind interpolation.
The *van Albada* scheme#footnote[Gale van Albada, _et al. A comparative study of computational methods in cosmic gas dynamics_, 1982.] is described by the limiter function

$ T B D $
It uses much less downwind interpolation than the van Leer scheme since $T B D$ as $T B D$ . From the diagram, it appears like a smooth-functioned, less diffuse version of minmod.


== Limiting multiple components

The calculation and application of a limiter is introduced in Sec. 3.11 for advection of a single scalar property. The figure below examines when the advected property is a vector, _e.g._ $T B D$ in the momentum conservation Eq. (2.47).

It shows a 2D cut along the $T B D$ - $T B D$ plane through a selection of cells in a regular mesh with the velocity vector $T B D$ displayed for each cell. When $T B D$ is interpolated to faces along the $T B D$ -direction using a limited scheme, the simplest approach is to calculate a limiter for each vector component, _e.g._ $T B D$ , and interpolate the component with that limiter.

In the example above, the profiles of $T B D$ and $T B D$ in the $T B D$ -direction are quite different, so the calculated limiters will be different for $T B D$ and $T B D$ components at each face. The limiting will *not* be _invariant_ under a rotation of the co-ordinate axes, leading to a different solution depending on the initial orientation of the geometry with respect to the axes.

The limiting can be invariant using a single limiter calculated from the magnitude $T B D$ , which is applied to all components of $T B D$ . The strength of the limiting corresponds to an average across the components, which is usually insufficient for the component which requires strongest limiting. This can cause instability.

Instead, the *‘V’ scheme* calculates the limiter based on the ‘worst-case’ direction, _i.e._ the direction of steepest gradient in $T B D$ at the cell face. It uses the following expression for $T B D$ for a vector $T B D$ , replacing Eq. (3.12) for a scalar:

$ T B D $

While the V scheme ensures invariance, it also provides greater stability than component limiting. It can remove oscillations in solutions, _e.g._ in the example above of supersonic flow over a step showing the effect in velocity in the cells adjacent to the step corner.

=== Multivariate limiting

Multivariate limiting applies the same limiter for advection discretisation across a set of 2 or more equations. It works by calculating the limiter $T B D$ for each solution variable in the equation set and applying the lowest $T B D$ to all equations.

It can be used in order to maintain consistency in the transport of individual fluid species, such as $T B D$ , $T B D$ , _e.g._ in the propagation of laminar flame (which is beyond the scope of this book).



== Linear upwind scheme

*Linear upwind* is another significant scheme for advection. It is used particularly for advection of momentum, _i.e._ $T B D$ or $T B D$ but can be effective for advection of other variables.

The scheme reduces the diffusive nature of upwind by including additional upwind cell values indirectly from the gradient $T B D$ in the upwind cell.

Linear upwind describes the face value as an extrapolation of the upwind cell value to the face using the upwind cell gradient $T B D$ and a vector $T B D$ from the cell centre to face centre. It first provides a contribution to the coefficients $T B D$ of a matrix equation $T B D$ by representing face values $T B D$ by the upwind value $T B D$ . The extrapolation is then introduced through an additional _explicit_ contribution $T B D$ to $T B D$ .

=== Skewness

Interpolation of values between cell centres is along a line joining the cell centres. Any interpolated value at a face relates to the point on the face intersected by that line.

_Skewness_ is the distance between the intersection point and the face centre. It can be represented by a vector $T B D$ or as the ratio $T B D$ , where $T B D$ is the distance from the face centre to edge in the direction of $T B D$ .

An interpolated value represents an average across a face, but as skewness increases, this representation becomes less accurate. High skewness (_e.g._ $T B D$ 1) does not immediately equate to poor accuracy, however, since interpolated values $T B D$ are multiplied by face areas within a surface integral, _e.g._ $T B D$ — and high skewness occurs at small faces.

Advection schemes do not generally include a correction for high skewness to improve accuracy. However, the linear upwind scheme naturally includes skewness correction since the explicit contribution is in the direction $T B D$ towards the face centre.



== Gradient discretisation

The discretisation of a gradient $T B D$ is exclusively an explicit calculation using current values of $T B D$ . The conservative form of gradient calculation is based on a surface integral.

From the gradient definition in Sec. 2.23, the discretisation is

$ T B D $
The face value $T B D$ is generally interpolated from cell values using the linear scheme.
=== Point linear interpolation

While skewness is generally not a concern for advection discretisation, it deserves greater attention in gradient calculation. For “bad” meshes, _e.g._ containing elongated tetrahedral cells, the point linear scheme is often adopted to reduce skewness error.

*Point linear* interpolation uses: the value $T B D$ , calculated using linear interpolation, which corresponds to the “face point” at the intersection of the line connecting cell centres and the face; and, values $T B D$ at each vertex, interpolated from adjacent cells using inverse-distance weighting.

The scheme breaks the polygonal face into triangles and calculates the average value $T B D$ at the 2 vertices and face point for each triangle, area $T B D$ . Point interpolation calculates the face value as the area-weighted average of triangle values, _i.e._ $T B D$ .

=== Least squares gradient

A gradient calculation using a *least squares* finite difference method is sometimes used within the finite volume framework. The method calculates the gradient in a cell which, when used to extrapolate the cell value to centres of all neighbouring cells, minimises the error between extrapolated values and cell values.

For a given cell, a tensor $T B D$ is calculated by summing over faces using the inverse distance weighting $T B D$ , where $T B D$ is the cell centre-centre vector:

$ T B D $
The gradient is then evaluated using the inverse of the tensor $T B D$ and values in the neighbour (N) and current (P) cells:
$ T B D $



== Gradient limiting

The standard gradient calculation can generate values with an unphysically large magnitude when cell quality is poor. A large gradient value in one cell can contribute sufficiently to the source term of an equation to cause unboundedness in the solution.

The calculated gradient can be _limited_ to avoid unphysical values. The limiting extrapolates the cell values $T B D$ using the gradient to _all_ neighbouring face centres. If the extrapolated $T B D$ at a face falls outside of the bounds of $T B D$ at neighbouring cell centres, the gradient is adjusted to match the bounding value.

This numerical scheme has an advantage that is takes into account all surrounding cells so is fully “multi-dimensional”. It ensures boundedness locally between cells so can help eliminate instabilities in _advection_ discretisation.

=== Linear upwind with gradient limiting

Limiting can be applied to the cell gradient in the explicit contribution of the linear upwind advection scheme, described in Sec. 3.14. The scheme with gradient limiting is analysed below for its TVD behaviour on a regular mesh in 1D.

In the following figure, the cell gradient is denoted by “ $T B D$ ” and the face gradient is denoted by “ $T B D$ ”. The ratio of gradients is $T B D$ from Eq. (3.12), which is adjusted by modifying $T B D$ . No limiting is applied for $T B D$ , when $T B D$ can be expressed by $T B D$ with $T B D$ . In that case, $T B D$ , corresponding to a function $T B D$ .

Limiting is applied for $T B D$ , reducing the gradient linearly to zero as $T B D$ . This corresponds to a function $T B D$ for $T B D$ . Similarly, limiting is applied when $T B D$ . Within that range, $T B D$ so $T B D$ .

On a Sweby diagram linear upwind with gradient limiting is TVD, with the line $T B D$ being the regime with no gradient limiting.

This function uses too much downwind interpolation to be a reliable TVD scheme for advection. However, as a gradient limiter for the linear upwind scheme, it is strongly bounded, due to the multi-dimensional nature of gradient limiting — with good accuracy due to the skewness correction in linear upwind.



== Time discretisation

A local time derivative $T B D$ in an equation can be discretised as a finite difference in time. Time $T B D$ is expressed in discrete intervals, or _steps_, of duration $T B D$ .

The *Euler* scheme calculates the derivative from the field at the current time $T B D$ and the previous, or _old_, time level $T B D$ by:

$ T B D $


Time discretisation thereby contributes to the _diagonal_ (_only_!) coefficients of the matrix $T B D$ and to the source vector $T B D$ .
=== Courant number

For a 1D domain in the $T B D$ -direction, the Courant number is the following dimensionless parameter for each cell of length $T B D$ :

$ T B D $
The Courant number originates from the solution of 1D advection, _e.g._ Eq. (2.32), with the Euler time scheme of Eq. (3.21) and an _explicit_ upwind advection scheme from Sec. 3.10. For that case, the Courant-Friedrichs-Lewy (CFL) condition#footnote[after Richard Courant, Kurt Friedrichs and Hans Lewy, _Uber die partiellen Differenzengleichungen der_ _mathematischen Physik_, 1928.] for convergence is $T B D$ across _all cells_. $T B D$ corresponds to a fluid particle moving across one cell in one time step, so its relevance to solution convergence is perhaps unsurprising.
In an _explicit_ solution, the convergence limit can reduce further to $T B D$ or $T B D$ , with more accurate schemes for advection. But the finite volume method is generally _implicit_ so stability can then be maintained with a higher maximum value of $T B D$ . Temporal accuracy is then the important consideration when choosing the $T B D$ of a simulation, both in terms of its mean and maximum value across all cells in the domain.

It is therefore important to monitor Courant number which needs to be calculated for 3D problems. Explicit, Euler, upwind discretisation of advection can be presented in 3D as follows:

$ T B D $
Here $T B D$ are positive fluxes which transports $T B D$ out of the cell of interest and $T B D$ are negative fluxes which transports $T B D$ from neighbouring cells. Since the CFL condition, $T B D$ , requires that the coefficient of $T B D$ cannot be negative, it follows that
$ T B D $
This 3D $T B D$ represents the volume of fluid leaving the cell in one time step, as a fraction of the cell volume, as shown above.


== Second order time schemes

A Taylor’s series expansion between the current time $T B D$ and ‘old’ time at $T B D$ relates to the Euler _implicit_ scheme by

$ T B D $
The truncation error $T B D$ , _i.e._ it is _first order accurate_ in time when the time derivative relates to $T B D$ _at the current time_, as shown in the figure below. Despite its low order, the Euler scheme is sufficiently accurate for many simulations since $T B D$ is generally small when it corresponds to $T B D$ .
Nevertheless, a second order time scheme may be required for simulations that demand higher temporal accuracy or to enable greater computational efficiency by running with larger $T B D$ .

=== Backward scheme

In Eq. (3.25) we can replace $T B D$ by the values $T B D$ at ‘old-old’ time $T B D$ . Subtracting the expression from Eq. (3.25) and rearranging terms gives the following relation for the _second_ derivative

$ T B D $
Substituting Eq. (3.26) into Eq. (3.25) gives the *backward* scheme which is second order accurate, using values at three time levels $T B D$ , $T B D$ and $T B D$ :
$ T B D $

=== Crank-Nicolson scheme

An implicit solution expresses the terms in an equation, _e.g._ advection, Laplacian, at the _current_ time. The *Crank-Nicolson* method,#footnote[John Crank and Phyllis Nicolson, _A practical method for numerical evaluation of solutions of partial_ _differential equations of the heat conduction type_, 1947.] expresses the terms at the _midpoint_ between the current and old times, to make the Euler time scheme second order accurate. Denoting discretised terms except the time derivative by $T B D$ , the Crank-Nicolson method solves

$ T B D $
where $T B D$ is _calculated_ using old time values $T B D$ .
A modern version of the scheme replaces the two $T B D$ factors by $T B D$ and $T B D$ , introducing the ‘offset coefficient’ $T B D$ , where $T B D$ corresponds to Euler implicit and $T B D$ is Crank-Nicolson Eq. (3.28). If $T B D$ is discretised _implicitly_ (as normal), the Crank-Nicolson scheme can be represented as a time derivative discretised by

$ T B D $
In practice, $T B D$ is generally used to ensure solution stability.


== Calculated derivatives

The introduction to matrix construction, Sec. 3.5, described implicit and explicit discretisation of terms in an equation. It concluded that the principal derivatives of $T B D$ that can be treated implicitly — forming matrix coefficients in $T B D$ — are a time derivative, advection and Laplacian.

Terms with other derivatives must be calculated from respective fields, _e.g._ $T B D$ from current values of $T B D$ . In Sec. 3.15, we have described discretisation of a gradient, which is always explicit. This section gathers together the other derivatives found in equations for fluid dynamics and associated models.

=== General divergence term

A general divergence term is any term that can be represented by $T B D$ . It excludes the Laplacian term which includes a gradient $T B D$ , and advection which includes $T B D$ .

The discretisation of a general divergence term is an explicit calculation using current values of $T B D$ . It is based on a surface integral using the divergence definition in Sec. 2.23 as shown below:

$ T B D $
The face value $T B D$ is generally interpolated from cell values using the linear scheme. Terms discretised using this scheme include $T B D$ in Eq. (2.45), a divergence of stress $T B D$ , _etc._
=== Curl of a vector

The curl derivative $T B D$ is calculated from the gradient $T B D$ and applying the Hodge dual operator given by Eq. (2.40) using the following relation:

$ T B D $
In other words, $T B D$ is discretised according to a scheme from Sec. 3.15, from which $T B D$ is calculated by Eq. (3.31).
=== Mag-square grad-grad

A derivative which appears in some model equations is $T B D$ , described as “mag-square grad-grad”. This derivative returns a scalar since the mag-square, _e.g._ $T B D$ , represents the inner product of $T B D$ with itself, as shown in Eq. (2.7).

The mag-square calculation always uses the appropriate inner product to reduce the result to a scalar. For a tensor $T B D$ , it is the double inner product, _i.e._ $T B D$ .

The grad-grad operator $T B D$ yields a third-rank tensor in the case that $T B D$ is a vector field. To avoid storing third-rank tensors, the mag-square grad-grad operator is evaluated by summing the result from the operator on each _component_ $T B D$ of $T B D$ by

$ T B D $
where $T B D$ is the number of components in $T B D$ .


== Other terms

In an equation for $T B D$ , there can be terms other than the time derivative, advection and Laplacian of $T B D$ , including:

- a linear function $T B D$ where $T B D$ is a scalar coefficient or field;
- a term without $T B D$ , sometimes including the derivative of another variable, _e.g._ $T B D$ .


The second example is simply discretised as an _explicit_ gradient term as described in Sec. 3.15. Derivative terms like this one are described in earlier sections and require no further discussion.

The first example, the $T B D$ term, requires much more discussion, particularly relating to the possibility of _implicit_ discretisation. Let us consider the following equation:

$ T B D $


If $T B D$ is discretised implicitly, then the matrix would contain a zero diagonal coefficient when $T B D$ , making it _singular_ or _non-invertible_. If this occurs, $T B D$ is not present in the linear equation for the relevant cell, so it cannot be solved.
Therefore, such a term _must_ be discretised _explicitly_ when it has a _negative sign_ (or positive on the right side of “=”), to ensure the matrix equation is solvable. The nature of Eq. (3.33) is that $T B D$ can only increase from an initial positive value.

=== Implicit discretisation of linear terms

Let us now consider the equivalent equation with a linear term that has a positive sign, _i.e._

$ T B D $
With Eq. (3.34), $T B D$ can only decrease from an initial positive value but reaches a lower limit at $T B D$ . Like many scalar properties, _e.g._ $T B D$ , $T B D$ , $T B D$ , _etc._, $T B D$ may have a lower physical bound of 0, which the equation intentionally reflects.
It is important to maintain a physical bound of $T B D$ . Discretisation of Eq. (3.34) using the Euler time scheme Eq. (3.21) gives

$ T B D $
which ensures boundedness since $T B D$ . An equivalent explicit discretisation gives $T B D$ which is only bounded when $T B D$ , similar to the $T B D$ limit imposed by explicit discretisation of advection, described in Sec. 3.17. To avoid this $T B D$ limit, we apply _implicit_ discretisation to terms with a _positive sign_ on the left side of “=”, where possible.
Even when the term is not linear in $T B D$ , it can be implemented as such by “dividing and multiplying by $T B D$ ”. For example, the $T B D$ turbulence model includes an $T B D$ term in the $T B D$ equation, _i.e._ $T B D$ , ignoring other terms. Dividing and multiplying the $T B D$ term by $T B D$ gives

$ T B D $
producing a linear term in $T B D$ which can be discretised implicitly using a coefficient $T B D$ .


== Terms which change sign

An equation may include a function $T B D$ which can return a value which is positive or negative, _e.g._

$ T B D $
Following Sec. 3.20, we aim to avoid a singular matrix in the case $T B D$ and ensure the solution exceeds the lower bound of 0 when $T B D$ . The discretisation of $T B D$ is therefore treated implicitly or explicitly within each cell based on
$ T B D $

=== Maintaining lower and upper bounds

The solution variable $T B D$ could have a physical bound at any low value $T B D$ and/or high value $T B D$ . We can treat the function $T B D$ in Eq. (3.37) as follows to maintain boundedness:

$ T B D $
Here, $T B D$ and $T B D$ are calculated using the current values of $T B D$ .
If we examine the situation with $T B D$ , Eq. (3.37) would be discretised as

$ T B D $
When we discretise $T B D$ implicitly alongside the Euler time scheme, Eq. (3.40) becomes:
$ T B D $
The solution for $T B D$ decreases when $T B D$ , but cannot fall below the $T B D$ bound within a solution step. In the limit that $T B D$ , the solution for $T B D$ no longer decreases.
An equivalent analysis shows Eq. (3.40) maintains boundedness for the upper bound $T B D$ . In that case, note that $T B D$ is negative so the term in $T B D$ can still be treated implicitly.

=== Fields bounded by 0 and 1

Many properties are expressed as a fraction, _e.g._ concentration of a chemical species, so are bounded by 0 and 1. When that happens, Eq. (3.39) simplifies to

$ T B D $
Here, $T B D$ and $T B D$ are calculated using the current values of $T B D$ .
It is particularly important to maintain boundedness of a field with 0-1 bounds when it is used to calculate another property which has a physical bound. For example, in a system of $T B D$ fluid phases, fluid $T B D$ can be calculated by $T B D$ from phase fractions $T B D$ and densities $T B D$ .

A small amount of unboundedness in $T B D$ , in the phase with highest $T B D$ , causes large unboundedness in $T B D$ , _e.g._ in 2 phases with $T B D$ and $T B D$ , the calculated $T B D$ with $T B D$ .



== Bounded advection discretisation

Sec. 2.9 described the conservation and boundedness qualities of the advective derivatives $T B D$ and $T B D$ , respectively. According to Eq. (2.31), the terms are equivalent when $T B D$ , _i.e._ when mass conservation Eq. (2.8) is satisfied in the case when $T B D$  constant.

During computation of an incompressible flow, numerical error can be significant, in particularly at the beginning of a steady-state calculation. The governing equations are not satisfied, including mass conservation, so $T B D$ .

Discretisation of an advection term of the form $T B D$ as described in Sec. 3.9 can cause unboundedness in $T B D$ . If $T B D$ has a physical bound that is violated, the solution rapidly becomes unstable.

The addition of the $T B D$ term ensures boundedness while $T B D$ , at the expense of conservation. Once $T B D$ , the solution is both bounded and conservative.

$ T B D $


This additional term is implicit in $T B D$ , so contributes $T B D$ to the diagonal coefficients of the matrix. If $T B D$ is positive, we might expect a singular matrix, as discussed in Sec. 3.20. However, the decrease to the diagonal coefficient is matched by an increase from the discretisation of $T B D$ . If that term uses the upwind scheme, the discretisation can be split into a contribution from positive outgoing fluxes $T B D$ and negative incoming fluxes $T B D$ . In extensive form (_i.e._ scaled by $T B D$ , see Sec. 3.24), using values in the cell of interest $T B D$ and neighbouring cells $T B D$ , the terms are
$ T B D $
The $T B D$ terms cancel, leaving negative coefficients for neighbouring cells (since $T B D$ is negative). The diagonal coefficient is then equal to the sum of magnitude of those neighbour cell coefficients, resulting in an invertible matrix.
=== Capturing physics

Boundedness and conservation can be compromised when a term in an equation does not capture the physics of the problem. An example from flame speed combustion modelling uses a parameter $T B D$ which represents the fraction of the unburnt fuel mixture.

The equation for $T B D$ includes the source term $T B D$ , where $T B D$ is a calculated flame speed. Its inclusion could cause the solution of $T B D$ to fall below its lower bound of 0.

Multiplying and dividing by $T B D$ changes the term to $T B D$ , where $T B D$ and unit vector $T B D$ . In this form, the term represents the non-conservative advection of $T B D$ by a flame velocity $T B D$ moving in the direction $T B D$ , which captures the physical nature of combustion. Boundedness can then be maintained by a suitable choice of advection scheme.



== Recommended discretisation schemes

The following figures provide a summary and recommendations of discretisation schemes for different equations. We start with momentum conservation Eq. (2.49), described in Sec. 2.13 and Sec. 2.14, for which _accuracy_ is the main consideration in scheme selection.



=== Equation for a bounded property

The discretisation of an equation for a bounded scalar field is illustrated using the $T B D$ equation from the $T B D$ turbulence model (assuming constant $T B D$ ). Schemes are selected with an emphasis on _boundedness_.





== Example of building a matrix equation

The previous sections describe methods to discretise derivative and other terms in order to build a matrix equation for a given physical equation. Let us demonstrate the construction of a matrix equation, using the momentum conservation equation from Sec. 3.23 as an example. It is a vector equation, so produces 3 matrix equations for $T B D$ , $T B D$ and $T B D$ .

The first term, the time derivative $T B D$ , might be discretised with the Euler scheme Eq. (3.21). Matrix equations are constructed in _extensive form_ as discussed in Sec. 3.6. Hence, the contributions from Eq. (3.21) to matrix coefficients $T B D$ and source vector $T B D$ are scaled by cell volume $T B D$ , _i.e._ $T B D$ and $T B D$ , respectively, as illustrated below.


The second term, the advective derivative $T B D$ , is discretised by Eq. (3.8). It pre-calculates the volumetric flux $T B D$ , using $T B D$ interpolated by Eq. (3.3) with linear weights Eq. (3.4).
The transported $T B D$ might be discretised using the linear upwind scheme described in Sec. 3.14. The scheme first applies upwind discretisation, which contributes outgoing positive fluxes $T B D$ to diagonal coefficients and negative fluxes $T B D$ to off-diagonals. It then adds an explicit contribution based on an extrapolated gradient $T B D$ (see Sec. 3.14). The gradient $T B D$ is usually calculated by Eq. (3.18) with gradient limiting from Sec. 3.16.


The third term, the Laplacian derivative $T B D$ , is discretised by Eq. (3.2). It requires $T B D$ , which is linearly interpolated from the cell centres. If the surface normal gradient $T B D$ includes a non-orthogonal correction $T B D$ , see Sec. 3.8, then the term contributes to $T B D$ and $T B D$ , as shown below.

The final term, $T B D$ is calculated using Eq. (3.18). Like all the other terms described here, it is implemented in extensive form, scaling by $T B D$ , so is calculated for each cell by the vector $T B D$ . The relevant component ( $T B D$ , $T B D$ , $T B D$ ) of this vector is then applied to the respective equation for $T B D$ , $T B D$ and $T B D$ .

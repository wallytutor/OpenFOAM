The previous chapters describe the components needed to build equations to solve problems in fluid dynamics, which are listed below:

- governing equations and common models, see Chapter 2;
- numerical methods to create discrete sets of linear equations, as matrix equations, see Chapter 3;
- boundary condtitions, which are applied to the matrix equations, see Chapter 4.


The final component is to assemble the set of matrix equations and apply _algorithms_ to:

- solve the set of linear equations described by each individual matrix equation;
- combine, or _couple_, solutions from multiple matrix equations to produce the overall solution.




This chapter first addresses the solution of individual matrix equations. It presents the most effective algorithms for solving matrix equations encountered with the finite volume method.

The important requirements of these _matrix solvers_ are robustness and efficiency, particularly for large simulations with many millions of cells. Solvers than accompany the finite volume method are generally _iterative_ in nature, in which the solution converges towards a prescribed acceptable level of accuracy, or _tolerance_, over a number of iterations.

The principal algorithms that couple sets of matrix equations will then be presented. The coupling of mass and momentum conservation with equations for $T B D$ and $T B D$ is challenging, and so requires particular attention.



== Structure of matrices

As described in Sec. 3.4, a matrix equation for solution variable $T B D$ contains of a set of coefficients $T B D$ where each row $T B D$ corresponds to the linear equation for the cell with index $T B D$ as follows:

$ T B D $

=== Matrix sparsity

Imagine creating a matrix equation for a transport equation for a scalar field, _e.g._ Eq. (2.65) with zero heat source $T B D$ .

The figure shows a matrix of size $T B D$ , where $T B D$ is the number of cells. Circles denote non-zero coefficients, which are filled for the diagonal coefficients ( $T B D$ ). The matrix $T B D$ is _sparse_, _i.e._ the majority of the coefficients $T B D$ are 0 (zero).

The sparsity is due to each cell interacting only with adjacent cells connected through its faces. For example, a 3D mesh of hexahedral cells produces up to 7 coefficients per matrix row, with one diagonal coefficient corresponding to a particular cell and 6 off-diagonal coefficients for the neighbour cells.

=== Matrix (a)symmetry

A symmetric matrix possesses the same coefficients across the diagonal, _i.e._ $T B D$ . The discretisation of a Laplacian term, _e.g._ $T B D$ , described in Sec. 3.7, produces coefficients that are symmetric because the surface normal gradient Eq. (3.5) uses the current and neighbour cell values in equal measure.

However, the discretisation of an _advection_ term, _e.g._ $T B D$ , generally produces _asymmetric_ coefficients. For example, if the upwind scheme is applied for flow from cell 1 into cell 2, then there would a contribution to $T B D$ but not to $T B D$ .

=== Matrix size

Parallel simulation allows affordable solution on huge meshes. For a mesh of size $T B D$ , the matrix would have $T B D$ coefficients, typically with $T B D$ that are are non-zero.

For reasons of efficiency, zero coefficients are not stored in the computer’s memory. Instead the storage provides an array of non-zero coefficients and addressing arrays of the corresponding row and column indices for each coefficient.



== Gauss-Seidel method

Finite volume numerics generally uses _iterative_ methods to solve each matrix equation. These methods calculate approximate solutions for $T B D$ , which become more accurate with successive iterative solutions.

Iterative methods are preferred because they are more efficient than direct methods, which solve a matrix equation exactly. Gaussian elimination, which is the numerical basis for direct solution methods, has a computational cost $T B D$ . This is prohibitive for many sizes of mesh in finite volume CFD.

*Gauss-Seidel*#footnote[Carl Friedrich Gauss, _letter to Christian Gerling_, 1823.] is a simple, iterative method which is generally effective for solving transport equations such as the example in Sec. 5.1. The method is illustrated by a sample equation

$ T B D $
The equations can be rewritten by: a) multiplying the off-diagonal coefficients by $T B D$ ; b) subtracting the result from the right hand side (r.h.s.); c) and, dividing by the diagonal coefficient. _i.e._:
$ T B D $

Starting with, $T B D$ , new values of $T B D$ are calculated by Eq. (5.3a), Eq. (5.3b) and Eq. (5.3c) in sequence, where the notation “ $T B D$ ” denotes “ $T B D$ is assigned the value of $T B D$ ”.

The first solution of Eq. (5.3a) is $T B D$ . The updated $T B D$ is substituted in Eq. (5.3b), whose solution is $T B D$ . Both updated values are substituted in Eq. (5.3c) to give $T B D$ .

The process is then repeated and through successive sweeps over the equations the solution converges as shown below.
Variable Start Sweep 1 …2 …3 …4 $T B D$ 0.0000 2.0000 2.4167 2.7431 2.8821 $T B D$ 0.0000 0.0000 0.5833 0.8069 0.9121 $T B D$ 0.0000 1.2500 1.6458 1.8392 1.9266 Variable Sweep 5 …6 …7 …8 …9 $T B D$ 2.9462 2.9755 2.9888 2.9949 2.9977 $T B D$ 0.9599 0.9817 0.9916 0.9962 0.9983 $T B D$ 1.9665 1.9847 1.9930 1.9968 1.9985
The error is $T B D$ , _i.e._ the difference between the approximate and exact values $T B D$ , $T B D$ and $T B D$ . After 9 sweeps $T B D$ for all variables, _i.e._ within  0.2% of the exact solution.

In summary, the Gauss-Seidel method is the following sequence of calculations for $T B D$ , repeated until convergence:

$ T B D $
The Gauss-Seidel method, _applied to sparse matrices_, has a computational cost $T B D$ , making it practical for simulations with large meshes that often occur in finite volume CFD.
Convergence of the method, and convergence measures for iterative methods in general, are discussed in the following sections.



== Convergence

The Gauss-Seidel method was demonstrated in Sec. 5.2 using a sample problem, Eq. (5.2), that converged to within 0.2% accuracy in 9 solution sweeps. This section discusses the criteria for convergence of a solution.

The easiest explanation of convergence considers the error $T B D$ for each equation $T B D$ , defined in Sec. 5.2. Substituting $T B D$ in each term in Eq. (5.3), _e.g._ in Eq. (5.3a), gives

$ T B D $
Since $T B D$ , Eq. (5.5) reduces to
$ T B D $
The magnitude of error $T B D$ is _at least_ as large as the sum of the terms on the r.h.s., but is smaller if the signs of $T B D$ and $T B D$ are different, _i.e._
$ T B D $
Repeating for Eq. (5.3b) and Eq. (5.3c) gives:
$ T B D $

The solution begins with an initial error of $T B D$ , $T B D$ and $T B D$ . After one sweep the error is $T B D$ , $T B D$ and $T B D$ .

The error is quickly distributed evenly, such that $T B D$ and $T B D$ are almost identical at sweep 2. The errors continue to reduce since Eq. (5.7) and Eq. (5.8) guarantee that no error is greater than the average of the other errors.

=== Condition for convergence

The behaviour of this problem indicates a convergence condition for the Gauss-Seidel method: the magnitude of the diagonal coefficient in each matrix row must be greater than or equal to the sum of the magnitudes of the other coefficients in the row; in one row at least, the “greater than” condition must hold.

This is known as _diagonal dominance_, which is a _sufficient_ condition for convergence, described mathematically as

$ T B D $
where the ‘ $T B D$ ’ condition must be satisfied for at least one $T B D$ . The description of the condition as “sufficient” means that convergence may occur when the condition is not met.


== Residual

In Sec. 5.3, we established a criterion for convergence of the Gauss-Seidel method. We now need a way to estimate a _level_ of convergence to determine when to stop iterating.

The analysis of convergence centred on the solution error $T B D$ , introduced in Sec. 5.2. In practice, $T B D$ cannot be determined since the exact solution is unknown. Instead the _residual_ provides a measure of the accuracy of the solution. The _residual vector_ $T B D$ represents the change to the solution of the equation, required to make $T B D$ exact, according to

$ T B D $
where $T B D$ . For the *remainder of this chapter*, the matrix notation $T B D$ is replaced by $T B D$ , equivalent to vector notation with an $T B D$ tensor $T B D$ .
The vector $T B D$ (of size $T B D$ ) provides _one value per matrix row_, with both positive and negative values. A measure of residual given by a _single_ value, is defined as

$ T B D $
where $T B D$ is the matrix norm, calculated as the sum of the magnitude of each component, _e.g._ $T B D$ ; the mean value of $T B D$ over all cells is denoted by $T B D$ .
The residual $T B D$ provides a measure of error in the _solution_ of $T B D$ , rather than the _absolute error_ $T B D$ . It is divided by the norms of $T B D$ and $T B D$ to reduce its dependency on the scale of the geometry and solution variable. By reducing its scale-dependency, $T B D$ can be used to compare the level of error equitably between simulations at different scales.

The figure above shows $T B D$ calculated from Eq. (5.11) following successive sweeps of the Gauss-Seidel method (starting from the initial $T B D$ ). The graph uses a logarithmic vertical scale since the values of $T B D$ extend over 4 orders of magnitude.

=== Tolerance

CFD software generally provides the following controls to stop the iterative solver:

- _absolute_ tolerance $T B D$ ;
- _relative_ tolerance $T B D$ .


Sweeping ceases if _either_ tolerance condition is satisfied: $T B D$ ; or $T B D$ , where $T B D$ is the initial residual within the particular solution step. The $T B D$ criterion is often deactivated by setting $T B D$ , especially for transient simulations when sufficient accuracy is required at every solution step.



== Diagonal dominance

Let us return to the transport equation in Sec. 5.1

$ T B D $

If the Gauss-Seidel method is applied to solve the equation, a convergent solution is only _guaranteed_ if the diagonal dominance condition of Sec. 5.3 is satisfied.

The contributions to the diagonal and off-diagonal coefficients are summarised below for each term in the equation.

The Laplacian term $T B D$ (note the sign), is discretised by Eq. (3.2) and Eq. (3.5), providing a positive contribution to the diagonal coefficient (relating to $T B D$ ) equal to the sum of the magnitude of negative off-diagonal coefficients (for $T B D$ ).

The contribution from the advection term $T B D$ , described by Eq. (3.8), is demonstrated above using the linear interpolation scheme of Eq. (3.4). The sign of $T B D$ depends on whether the flux in incoming or outgoing from the cell. Some contributions to off-diagonal coefficients are thus negative and some positive, with the diagonal contribution tending to zero as the mesh becomes more regular and flow becomes more uniform.

Advection with linear interpolation maintains diagonal equality while its positive off-diagonal contributions offset the negative contributions from the Laplacian term.

But once the advection contributions exceed the magnitude of the Laplacian contributions, diagonal equality is *not* achieved. This occurs when $T B D$ (= 2 for a regular mesh), where the “numerical” Péclet number $T B D$ is specified at cell faces from its definition in Sec. 2.21 with $T B D$ .

Diagonal equality is achieved with _upwind interpolation_ and any scheme that provides an explicit correction on upwind, _e.g._ linear upwind. This is conditional on mass conservation being satisfied, as discussed in Sec. 3.22. If mass is not conserved, the bounded advection scheme, described in that section, is required to ensure diagonal equality.

The time derivative $T B D$ increases the diagonal coefficient only, so promotes diagonal dominance. Notably, the diagonal contribution of $T B D$ from $T B D$ becomes larger than $T B D$ from advection with upwind when $T B D$ .

In summary, _diagonal dominance is not guaranteed_ due to the contribution of the advection term.



== Under-relaxation

Sec. 5.5 concludes that the matrix of a typical transport equation is not guaranteed to be diagonally dominant. Some action may therefore be required to ensure a convergent solution.

Under-relaxation is a general method used to improve solution convergence by limiting the amount a variable changes during a solution step.

During a solution step, assume a single value of a field $T B D$ in one cell changes from its current value $T B D$ to the new value $T B D$ . Under-relaxation would limit the change $T B D$ by a fraction $T B D$ , $T B D$ , so that the value taken from that solution step is

$ T B D $
In some situations, Eq. (5.13) is applied after a solution step. This simple approach is known as _field under-relaxation_, which has one notable disadvantage that it requires additional storage of the intermediate field $T B D$ in computer memory.
When a solution step involves solving a matrix equation, the new values $T B D$ come from an iterative method like Gauss-Seidel. Combining the under-relaxation of Eq. (5.13) with the Gauss-Seidel calculation of Eq. (5.4) gives:

$ T B D $
Rearranging Eq. (5.14) gives the following relation:
$ T B D $
Equation 5.15 is simply the matrix equation $T B D$ modified by:
- increasing the diagonal coefficients $T B D$ by division by $T B D$ ;
- multiplying the difference between the new and original $T B D$ coefficients by the current $T B D$ and adding it to the source $T B D$ .


Modifying the matrix equation this way, known as _equation under-relaxation_, provides an alternative to Eq. (5.13) for under-relaxing a solution of $T B D$ , without the temporary storage of $T B D$ .

=== Ensuring diagonal dominance

The modification to $T B D$ expressed by Eq. (5.15) inspires a strategy to ensure diagonal dominance of the matrix as follows.

Each diagonal coefficient which does not satisfy Eq. (5.9) is increased until it is diagonally equal. The change to the coefficient is multiplied by the current $T B D$ and added to $T B D$ .

This approach to ensure diagonal dominance is effective since it only modifies matrix coefficients where necessary. Otherwise, if the discretisation schemes, and $T B D$ and $T B D$ are favourable, then no changes to the matrix are necessary.



== Iterative solution

Imagine a problem with steady flow in which $T B D$ is changing at a boundary, _e.g._ at an inlet. The system could be simulated by solving Eq. (5.12) with a $T B D$ field which is constant in time.

The simulation would evolve by an iterative sequence which increases time $T B D$ by an increment $T B D$ , then solves the equation, as shown below.

=== Setting tolerances

The $T B D$ equation is solved using an iterative method such as Gauss-Seidel. Tolerance parameters determine when the method stops iterating within the current solution step.

To capture the transient solution for $T B D$ accurately, the equation would need to converge to a suitable absolute tolerance $T B D$ at every time step. The question is: what $T B D$ is suitable to give an accurate solution? — without setting it to a very small value, _e.g._ $T B D$ , which would would be highly inefficient, due to an excessive number of solver iterations.

Solution accuracy should be determined by some measure, or _metric_, _e.g._ temperature rise, pressure drop, drag coefficient, _etc._, relevant to the _objective_ of the simulation.

The metric should be initially calculated using a value of $T B D$ , _e.g._  $T B D$ . It can then be calculated again at a lower $T B D$ , typically decreasing by one order of magnitude, _i.e._ $T B D$ .

The metrics are then compared. If the difference is within the required level of accuracy, then $T B D$ is optimal. Otherwise, $T B D$ would need to be decreased further and the process repeated.

=== Residual calculation at intervals

The residual calculation of Eq. (5.11) incurs a significant computational cost since it involves multiplications of matrices. That cost can be reduced by setting an interval on the _number of_ _sweeps_ between residual calculations.

The downside is that convergence checks can only happen when the residual is calculated. If the residual is calculated at intervals of every second sweep, a solution that would otherwise converge at an odd number of sweeps is prolonged by a extra sweep.

Residuals are high when a simulation starts so a relatively high number of sweeps, _e.g._ 10, may be needed for convergence. The saving of calculating residuals at an interval may then offset the additional cost of unnecessary additional sweeps, giving a lower net cost per time step. However, as the residuals decrease and fewer sweeps are needed, attempts to reduce cost by intermittent residual calculation can be ineffective or counter-productive.



== Accelerating convergence

The Gauss-Seidel method sweeps across cells, calculating the changes in values of fields, _e.g._ $T B D$ . Each sweep propagates the changes in $T B D$ through cells in the domain in much the same way that the flow transports $T B D$ physically.

The method is therefore remarkably effective, so the scope to accelerate convergence to reduce computational cost is limited. Two options that may reduce cost, relating to mesh numbering, are described in this section.

=== Symmetric Gauss-Seidel

The Gauss-Seidel method sweeps through equations in sequence $T B D$ of cell numbering. The changes in the solution “flow” in a direction, based on the choice of cell numbering which may be arbitrary and therefore not optimal.

The *symmetric Gauss-Seidel* method follows the standard $T B D$ sequence initially, known as a _forward_ sweep. For every second sweep, the sequence is executed in reverse, _i.e._ $T B D$ , known as a _backward_ sweep.

The change in sweep direction moderates the directional bias of the standard Gauss-Seidel method which may reduce the total iterations for convergence and reduce cost.

=== Mesh numbering

When a matrix equation is solved in sequence $T B D$ , values from one row of the equation are substituted in subsequent equations. Convergence improves if the substitutions are made quickly, _e.g._ in the subsequent 2 or 3 equations rather than in the last equation.

Convergence is therefore affected by the cell numbering which determines the position of non-zero coefficients in the matrix. If the numbering is approaching random, _e.g._ as a consequence of the mesh generation, convergence is sub-optimal.

Where this occurs, a method such as reverse Cuthill-McKee#footnote[Elizabeth Cuthill and James McKee, _Reducing the bandwidth of sparse symmetric matrices_, 1969; “reversed” by Alan George and Joseph Liu, _Computer Solution of Large Sparse Positive Definite Systems_, 1981.] can be applied to renumber the cells of a mesh to cluster non-zero coefficients around the diagonal, as shows above.

The matrix in Sec. 5.1 is from a mesh generated in a structured manner with ordered indexing along rows of cells. Renumbering a matrix like that generally provides little benefit since it already has many non-zero coefficients clustered around the diagonal.



== Systems of equations

Most CFD calculations involve solving a system of equations that represent the physics of the problem. For example, laminar flow by natural convection can be represented by the equations introduced in Sec. 2.20, reproduced below.

The system provides 3 equations (1 vector, 2 scalar) which can be solved for 3 unknowns, $T B D$ , $T B D$ and $T B D$ .

As discussed in Sec. 3.4, the finite volume method solves an individual matrix equation for each variable, _e.g._ $T B D$ for $T B D$ . The vector equation for $T B D$ is decoupled into 3 matrix equations for individual components, _i.e._ $T B D$ , $T B D$ and $T B D$ .

Each individual matrix equation for one solution variable, _e.g._ $T B D$ , incorporates terms from other variables, _e.g._ $T B D$ , into the source vector $T B D$ . The contribution to $T B D$ is calculated using current values of the respective variables. Systems of equations are thereby solved by _successive substitution_ of solved variables into the source vectors of subsequent equations.

An iterative solution for a single equation, like the one in Sec. 5.7, can be extended to a system of equations. Time $T B D$ is incremented by $T B D$ and equations are solved in sequence, before returning to start the next time step with the increment of $T B D$ .

The substitutions in the momentum and pressure equations are particularly important, culminating in _corrections_ to $T B D$ and the advective flux $T B D$ , discussed in Sec. 5.10.

At the start of any time step the current $T B D$ becomes $T B D$ for the discretisation of the $T B D$ term in the momentum equation. The advection term $T B D$ is discretised by Eq. (3.8), treating one $T B D$ as flux $T B D$ and the other as the advected quantity.

The equation is solved for $T B D$ . The new solution for $T B D$ is substituted into the $T B D$ equation which is solved for $T B D$ . The new solution for $T B D$ is then used to _correct_ $T B D$ in order to help enforce the mass conservation constraint $T B D$ ( $T B D$ ).

Before the current solution step is completed, $T B D$ is also corrected to reduce the error in the discretisation of $T B D$ when it then becomes $T B D$ in the following solution step. The correction also provides a better “initial guess” $T B D$ for the matrix solution of the next momentum equation, which reduces the solution time.



== Pressure-velocity coupling

The previous section combined equations for $T B D$ and $T B D$ , governing momentum and mass conservation, in a sequential solution.

The algorithms used to couple these equations, in a manner which is convergent, uses the following notation to describe terms in the momentum equation, _e.g._ Eq. (2.67), _excluding_ $T B D$

$ T B D $
where:
- $T B D$ is a linear term in $T B D$ ;
- $T B D$ is a function of $T B D$ and other sources.


=== Momentum corrector

A _momentum corrector equation_ is formed by expressing the momentum equation, _e.g._ Eq. (2.67) in terms of the notation of Eq. (5.16), and rearranging to give



$ T B D $
The equation provides an update to $T B D$ , based on current values of $T B D$ and $T B D$ substituted on the r.h.s. In other words, $T B D$ and $T B D$ are calculated explicitly.
In the algorithms described in the following sections, $T B D$ is chosen to be the extracted diagonal coefficients $T B D$ $T B D$ of the matrix $T B D$ corresponding to the momentum terms in Eq. (5.16).

The split between $T B D$ and $T B D$ is shown above. As a matrix with diagonal components only, $T B D$ has one value per cell so can be represented as a scalar field. Setting $T B D$ to be the diagonal coefficients is a natural choice for implicit treatment of $T B D$ within the coupling algorithm.

=== Flux corrector

A _flux corrector equation_ follows from Eq. (5.17) by interpolating $T B D$ to cell faces and evaluating $T B D$ according to

$ T B D $

=== Pressure equation

A _pressure equation_ is then created by substituting fluxes $T B D$ from Eq. (5.18) into the mass conservation Eq. (2.46) in discrete form $T B D$ . The resulting expression is equivalent to a discretised pressure equation, with coefficients containing $T B D$ and $T B D$ ,

$ T B D $



== Boundary fluxes

The equations in Sec. 5.10 are used to compute the velocity $T B D$ , flux $T B D$ and pressure $T B D$ based on conservation of mass and momentum. At boundaries, the flux corrector Eq. (5.18) must calculate $T B D$ in a manner consistent with $T B D$ and $T B D$ and their respective boundary conditions.

For example, at an impermeable stationary wall the calculated flux must be $T B D$ , consistent with the no-slip condition $T B D$ . At boundary faces, $T B D$ and $T B D$ must be compatible to evaluate the correct $T B D$ according to Eq. (5.18).

The figure shows the fundamental boundary conditions from Sec. 4.3 and corresponding flux evaluations. At the inlet and wall boundaries, $T B D$ is directly assigned from the boundary velocity $T B D$ by $T B D$ within the flux corrector Eq. (5.18).

In the absence of body forces, the $T B D$ boundary condition is commonly applied, as discussed in Sec. 4.4. The flux $T B D$ is then equivalent to assigning $T B D$ in Eq. (5.18).

At the outlet, $T B D$ is not prescribed since $T B D$ is not a fixed value condition. Instead, $T B D$ is evaluated from Eq. (5.18) using $T B D$ taken from cells adjacent to the boundary and $T B D$ calculated on the boundary.

=== Fluxes with a body force

When a body force $T B D$ is present in the momentum equation, the gradient condition for $T B D$ , _e.g._ at an inlet or wall, in principle becomes $T B D$ , as discussed in Sec. 4.4. The precise details of the boundary condition in fact depend on how $T B D$ is incorporated within the coupling algorithm, discussed below.

One approach, illustrated by the algorithm in Sec. 5.10, is to include the body force $T B D$ within $T B D$ in Eq. (5.16). In that case, the assignment $T B D$ cannot be valid, so $T B D$ is adopted instead.

With $T B D$ established, the $T B D$ condition is calculated based on the known $T B D$ by inverting Eq. (5.18). This approach causes $T B D$ to include a contribution from viscous stresses, as in Eq. (4.5), which may cause instability.

To avoid this problem, $T B D$ is omitted from $T B D$ in Eq. (5.16), appearing instead as an extra term in the other equations in Sec. 5.10, _e.g._ the flux corrector Eq. (5.18) which becomes

$ T B D $
Assigning $T B D$ at boundaries where $T B D$ is fixed satisfies Eq. (5.20) when $T B D$ . This is the condition ignoring viscous stresses described in Sec. 4.4.


== Steady-state solution

The equations in Sec. 5.10 are combined using algorithms that couple the solutions for $T B D$ and $T B D$ . One algorithm is SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)#footnote[Suhas Patankar and Brian Spalding, _A calculation procedure for heat, mass and momentum transfer in_ _three-dimensional parabolic flows_, 1972.], which is presented here with a modern interpretation.

The SIMPLE algorithm is generally used to produce steady flow solutions in CFD. These solutions are directly applicable for flows that reach a _steady state_, _i.e._ when flow variables stop changing in time. They can also provide approximate solutions to flows that are moderately unsteady, usually at a lower cost than a more exact transient solution.

An example of the algorithm is shown for the system of equations presented in Sec. 5.9. The time derivative ( $T B D$ ) terms are omitted due to the steady-state assumption.

The algorithm involves an iterative sequence with steps, $T B D$ . It begins by constructing a matrix equation for energy which is under-relaxed by a factor $T B D$ . The equation is solved for $T B D$ , which is used to update $T B D$ according to an equation of state. A matrix equation is then constructed using all the terms from the momentum equation _excluding_ $T B D$ , _i.e._

$ T B D $
The matrix equation is under-relaxed by a factor $T B D$ before equating with $T B D$ and solving for $T B D$ (the _momentum predictor_).
$T B D$ and $T B D$ are then evaluated from $T B D$ (the _momentum matrix_), as described on page 351. They are used to form the _pressure_ _equation_, which is solved for $T B D$ .

The new pressure $T B D$ is used to correct the flux $T B D$ so that it obeys mass conservation better (the _flux corrector_). It is then under-relaxed by a factor $T B D$ before correcting $T B D$ before the next solution step begins (the _momentum corrector_).



== Steady-state convergence

The absence of a time derivative from an equation written in steady-state form reduces the diagonal dominance, and hence convergence, of the resulting matrix equation, as discussed in Sec. 5.5. Under-relaxation is therefore applied to $T B D$ $T B D$ , and $T B D$ to promote convergence in the algorithm in Sec. 5.12.

The $T B D$ and $T B D$ fields use equation under-relaxation with factors $T B D$ and $T B D$ , respectively. A value of 0.7 is commonly applied, decreasing to 0.5 for less convergent cases (and sometimes to 0.3 in compressible flow cases, beyond the scope of this book).

Under-relaxation of $T B D$ is more subtle. The flux corrector requires that $T B D$ is not under-relaxed to ensure $T B D$ obeys mass conservation better. For the momentum corrector, field under-relaxation is subsequently applied to $T B D$ with a factor $T B D$ . To find an optimal $T B D$ for the momentum corrector, we examine

$ T B D $
This explicit momentum equation contains diagonal coefficients $T B D$ and an off-diagonal contribution $T B D$ from neighbour cell coefficients $T B D$ and associated velocities $T B D$ .
Convergence is compromised by the explicit nature of $T B D$ . It can be more implicit by “adding and subtracting” $T B D$ , where coefficients $T B D$ are applied to $T B D$ in the “owner” cells:

$ T B D $
Combining Eq. (5.22) and Eq. (5.23) gives
$ T B D $
which attempts to make Eq. (5.22) more implicit in $T B D$ .#footnote[Jeff Van Doormaal and George Raithby, _Enhancements of the SIMPLE method for predicting_ _incompressible fluid flows_, 1984.] The pressure equation derived from Eq. (5.24) is:
$ T B D $
This corresponds to under-relaxation of Eq. (5.19) by $T B D$ , where $T B D$ . A momentum matrix with coefficients $T B D$ is approximately diagonally equal since there is no time derivative. Since $T B D$ represents the diagonal coefficients under-relaxed by $T B D$ , $T B D$ . Relating the expressions for $T B D$ and $T B D$ gives the optimal under-relaxation factor for $T B D$ for convergence as
$ T B D $
leading to the popular choice of $T B D$ and $T B D$ .
=== Residual control

The algorithm in Sec. 5.12 uses a fixed number of solution steps $T B D$ . In practice, $T B D$ must be chosen to be large enough to reach an acceptable level of convergence. Once convergence is reached, the simulation should stop to avoid unnecessary computing cost.

A common stopping criterion applies a residual level for each equation, below which the equation is deemed to be converged. When all equations satisfy their respective residual controls, the simulation then stops.

Convergence can also be determined by monitoring any suitable metric, including objective measurements from the simulation, _e.g._ a force coefficient. When the metric no longer changes significantly over subsequent steps, the simulation is stopped.



== Descent methods

The Gauss-Seidel method, introduced in Sec. 5.2, provides convergent solutions for many problems in CFD. It is most effective when a modest reduction in residual is required, _e.g._ as part of a steady-state solution described in Sec. 5.12.

When the Gauss-Seidel method requires a lot of sweeps (_e.g._ over 10) to converge to a suitable tolerance, alternative methods may be more efficient. _Descent methods_ provide alternative matrix solvers that are often used in CFD.

Descent methods represent the equations which are being solved as a _minimisation_ problem. It is demonstrated below using a matrix equation of the form $T B D$ with 2 values

$ T B D $
The minimisation presents the equation in _quadratic form_; in matrix notation, it is a _scalar_ function of the form
$ T B D $
The quadratic form of Eq. (5.27), is:
$ T B D $
The gradient of $T B D$ is:
$ T B D $
Critically, $T B D$ , _i.e._ the negative of the residual vector Eq. (5.10), as verified by the model example.
Equating the gradient to zero, $T B D$ , corresponds to a _minimum in_ _the quadratic function_. At the same time, it is the _solution_ to $T B D$ . The method is therefore concerned with finding the minimum of the quadratic form efficiently.

For this method to work, the quadratic form must have a minimum, which requires that $T B D$ is _symmetric_ and _positive-definite_. A positive-definite matrix is hard to visualise, but for a 2-value function it ensures the quadratic function is a paraboloid.

Diagonal dominance is the convergence condition for the Gauss-Seidel method, discussed in Sec. 5.3. Importantly, a symmetric matrix that is diagonally dominant, and has positive diagonal coefficients, is also positive-definite.

=== Matrix operations

The ‘ $T B D$ ’ operation between two single-column matrices, _e.g._ $T B D$ , in Eq. (5.28) is represented in other texts using matrix notation by a transpose $T B D$ . For the example in Eq. (5.27), it is:

$ T B D $



== Conjugate gradient method

The solution to a matrix equation $T B D$ by descent methods involves finding the minimum of the equation in quadratic form. This can be illustrated using the contours of the paraboloid describing the quadratic for the case with values $T B D$ and $T B D$ .

The search for the minimum involves a series of updates to $T B D$ of the form

$ T B D $
where $T B D$ is the next update using latest (current) values $T B D$ . The column vector $T B D$ provides the _direction_ of the line search towards the minimum; the scalar $T B D$ provides the magnitude of the line is that direction.
=== Steepest descent

The intuitive way to reach the minimum is to follow the direction of _steepest descent_. The method is fairly simple because the direction is defined by the negative of gradient of the quadratic form $T B D$ which is the residual $T B D$ .

The distance to “walk” is naturally until the lowest point is reached, corresponding to#footnote[Jonathan Richard Shewchuk _An introduction to the conjugate gradient method without the agonizing_ _pain_, 1994.]

$ T B D $
The method reaches the minimum point along a zigzag path where consecutive directions are orthogonal. The directions do not change and result in a _larger_ number of very short steps, the _more_ the paraboloid is stretched in one direction.
=== The conjugate direction

The *conjugate gradient* (CG) method chooses search directions that are _conjugate_ with $T B D$ . This means each new direction $T B D$ corresponds to the previous one $T B D$ , satisfying

$ T B D $

This can be imagined as the directions being orthogonal _with the stretch in the paraboloid_ _removed_. For 2 values, CG finds the minimum in 2 steps, rather than several zigzag steps.

CG provides the basis for practical matrix solvers for CFD, described in Sec. 5.16. For a detailed explanation of the CG method, see the reference below.



== Preconditioning and asymmetry

The conjugate gradient (CG) method provides a good basis for a matrix solver in CFD. In theory, it provides an exact solution in $T B D$ steps for a problem with $T B D$ values. In practice, there is an accumulation of rounding errors from finite precision arithmetic, such that an exact solution is not attained. Problems are also too large for it to be feasible to run $T B D$ steps.

CG is instead used as an iterative method, aiming to get the solution residual within an acceptable tolerance. The rate of convergence is therefore critical, which can be shown to be slower for increasing _condition number_ of the matrix, $T B D$ , where $T B D$ and $T B D$ are its largest _eigenvalues_.

In the example with two values, the eigenvalues correspond to the gradients of the principal axes of the paraboloid described by the quadratic function. Lower $T B D$ corresponds to closer gradients, _i.e._ to a rounder “bowl” shape.

The best convergence corresponds to $T B D$ when $T B D$ ( $T B D$ is the identity matrix). _Preconditioning_ finds a matrix $T B D$ , similar to $T B D$ , which can be inverted and multiplied to the equation as follows:

$ T B D $
The aim of preconditioning is to improve convergence by choosing a preconditioner $T B D$ so that $T B D$ is closer to $T B D$ than $T B D$ . Effective methods for calculating the preconditioner $T B D$ are:
- diagonal-based incomplete Cholesky (DIC) for symmetric $T B D$ ;
- diagonal-based incomplete LU (DILU) for an asymmetric $T B D$ .


=== Biconjugate gradient stabilised method

The *preconditioned CG* (PCG) method minimises the quadratic form which requires that $T B D$ is symmetric. However, any equation containing an advection term produces an asymmetric matrix (see Sec. 5.1), when PCG is unsuitable.

The *preconditioned biconjugate gradient* (PBiCG) method is designed for asymmetric matrices. It augments the calculations of residuals and line search directions in the CG method with a second set of calculations relating to the transpose matrix $T B D$ . The computational cost is approximately 2 $T B D$ the cost of PCG.

The convergence behaviour of PBiCG is somewhat erratic, with large variations in the reduction of the residual between successive iterations. Sometimes the method breaks down and fails.

The *preconditioned biconjugate gradient stabilised* method (PBiCGStab) is a modification of PBiCG. It is exhibits much smoother convergence and robustness than PBiCG so it is far more preferable for CFD applications.

Recommended CG-based methods are summarised below.

- _symmetric_: *PCG* with *DIC* preconditioning.
- _asymmetric_: *PBiCGStab* with *DILU* preconditioning.


The implementation of preconditioning and CG-based methods can be found in: Barrett R. _et al._ #link("http://www.netlib.org/templates/templates.pdf")[_Templates for the solution of linear systems_], available at #link("http://www.netlib.org")[http://www.netlib.org].



== Multi-grid method

The methods for solving matrix equations described so far are Gauss-Seidel and CG. They are iterative, and require few iterations to calculate changes in the field, when the change at any given point is influenced only by points in the close vicinity.

This makes them efficient for solving equations whose range of influence is limited, _e.g._ by the flow speed through advection $T B D$ in a transport equation.

The pressure equation contains neither a local rate of change $T B D$ nor advection term (assuming $T B D$ ). A disturbance at any point influences the solution everywhere in the domain instantaneously, as discussed in Sec. 4.3.

To transfer changes across a domain, Gauss-Seidel might require as many sweeps as the number of cells across a domain, which would be impractical. Some method is therefore needed that transfers change across the domain more efficiently.

The multi-grid method#footnote[Achi Brandt, _Multi-level adaptive solutions to boundary-value problems_, 1977.] uses a coarse mesh to overcome the problem of slow transfer of change across the domain, reducing the error at “low-frequency”. It transfers that solution through a succession of finer meshes to be accurate at higher frequencies.

The multi-grid method works by first calculating the matrix $T B D$ and residual vector $T B D$ on the original (finest) mesh. Coarser meshes are then successively formed until the mesh reaches its coarsest level, containing as few as 10 cells.

The matrix $T B D$ and $T B D$ are formed at each level of coarsening. Different methods exist to calculate $T B D$ and $T B D$ , including a simple _agglomeration_ discussed in Sec. 5.18.


Equations are constructed in terms of the _correction_ $T B D$ , required to make $T B D$ exact, related to the residual $T B D$ by
$ T B D $
At the coarsest level, $T B D$ is solved precisely for $T B D$ , which can be performed efficiently using a direct solver since the matrix is small.
$T B D$ is then _injected_ into the finer level, providing the initial $T B D$ for a few sweeps of the Gauss-Seidel method at that level. Smoothing and injection are repeated up to the finest level, when the final $T B D$ is applied to $T B D$ to yield the solution.

It is generally more efficient to do more Gauss-Seidel sweeps at the coarser levels, when the cost if low due to a smaller matrix size. For example, 4 sweeps may be applied at the coarser levels, while 2 sweeps are applied at the finest.



== GAMG method

*GAMG* is an effective form of multi-grid method used in finite volume CFD. It combines:

- geometric (“G”) agglomeration to define the structure of the coarse meshes;
- the _algebraic_ multi-grid (AMG) method, where the matrix $T B D$ is constructed at a coarser level from coefficients at the finer level, rather than by geometric data from the coarse mesh.


=== Agglomeration

_Pairwise agglomeration_ forms coarser meshes by joining pairs of cells at each level of coarsening. In a sweep through the cells, an unpaired cell is paired with the (unpaired) neighbour that shares the face with largest area. This method generally maintains a low aspect ratio (see Sec. 8.2) in the resulting agglomerated cells.


The example shows an agglomeration of cell pairs 1-2, 3-4 and 5-6, which form cells 1, 2 and 3 respectively in the coarse mesh. If we solve the equation $T B D$ on a 2D mesh of square cells with $T B D$ , then Eq. (3.2) and Eq. (3.5) calculate diagonal coefficients $T B D$ and non-zero off-diagonal coefficients $T B D$ .
*Algebraic multi-grid* creates coefficients in the coarser mesh by summing coefficients $T B D$ and source $T B D$ from the finer mesh. The example produces one row of $T B D$ which is $T B D$ .

If the Laplacian for the coarse mesh were discretised directly, the coefficients in that row would be $T B D$ due to increasing $T B D$ . This discrepancy between agglomerated and calculated coefficients is repeated in subsequent agglomerations.

Multi-grid begins solving/smoothing for $T B D$ (with $T B D$ ) at the coarsest mesh. Coarse cell values of $T B D$ are then injected into corresponding cells in the next finest mesh.

With a Laplacian term, the coefficient discrepancy caused by algebraic agglomeration causes $T B D$ to be under-predicted. A correction to $T B D$ can be applied by scaling it by a factor

$ T B D $
This scaling is used when the equation being solved is dominated by a Laplacian term, _e.g._ a pressure equation, when summing coefficients between agglomerated cells is inaccurate. Scaling is not generally used for transport dominated by advection, since there is no error associated with summing coefficients from the advective derivative.
Equation 5.37 is derived by minimising the error with respect to $T B D$ in the equation $T B D$ . The minimisation is performed by setting $T B D$ for the quadratic form

$ T B D $



== Transient solution

A modern interpretation of the SIMPLE algorithm was described in Sec. 5.12 to couple steady solutions for $T B D$ and $T B D$ . An equivalent _transient_ solution follows an iterative sequence in which equations for $T B D$ and $T B D$ are solved over successive time intervals $T B D$ , between a start time $T B D$ and end time $T B D$ .

In a transient simulation, $T B D$ needs to be relatively small to maintain sufficient accuracy as the solution evolves over time $T B D$ . The equations do not then generally require under-relaxation to converge due to the increase of $T B D$ to the diagonal coefficients of the matrix from discretisation of the time derivative $T B D$ .

A transient simulation could follow the same sequence as the SIMPLE algorithm on page 359, but iterating over that sequence within each time step is costly. A more efficient algorithm solves the sequence only once per time step but adds an iterative loop which substitutes $T B D$ from the momentum corrector into the $T B D$ term of the pressure equation.

The pressure equation is then solved a second time, followed by a second momentum corrector, in the style of the PISO algorithm.#footnote[Pressure Implicit with Splitting of Operators, 1986] The addition of this “PISO loop” improves the accuracy of $T B D$ , $T B D$ and $T B D$ within each time step. The improved $T B D$ becomes $T B D$ in the next time iteration, which critically increases the accuracy of the time derivative $T B D$ in the momentum equation.

Without any iteration over the entire system of equations, the advection terms are discretised using $T B D$ from the previous time step. This “lagging” of $T B D$ does not compromise accuracy significantly, if $T B D$ is sufficiently small, e.g. corresponding to $T B D$ .

A further iteration of the PISO loop can be introduced to solve a third pressure equation and momentum corrector, in particular as part of an update to non-orthogonality described in Sec. 5.20. Increasing the iterations still further is generally not beneficial.



== Transient solution controls

The following additional controls are commonly used for the transient solution described in Sec. 5.19:

- an iterative loop around the pressure equation to update any non-orthogonal correction;
- an on/off switch for the momentum predictor.


=== Non-orthogonal corrector loop

The algorithm solves the pressure equation, Eq. (5.19), which includes the discretisation of the Laplacian term $T B D$ , according to Sec. 3.7. For a computational mesh with significant non-orthogonality, the corrected scheme Eq. (3.7) should be applied to maintain sufficient accuracy.

This correction to the pressure Laplacian is calculated using the current $T B D$ , which produces an explicit contribution to the matrix equation. The equation can be imagined with the explicit correction as an additional term $T B D$ .

To improve accuracy within each time step during a transient simulation, $T B D$ can be recalculated using the latest $T B D$ and the equation solved again.

The accuracy can be improved by further iterations, but this is not generally recommended due to the increased computational cost. With the pressure equation encountered twice as part of the PISO loop, and with one corrector on each occasion, the pressure matrix is solved a total of 4 times.

Solving the pressure matrix is costly, so any reduction in the number of solutions is welcome. Recognising that $T B D$ is updated as part of the PISO iterations anyway, an alternative solution strategy uses no non-orthogonal correction but adds an extra PISO loop, solving the pressure equation a total of 3 times.

=== Solving the momentum predictor

The transient algorithm on page 385 includes the momentum predictor that provides an initial solution for $T B D$ , which is used to evaluate $T B D$ in the pressure equation. The solution for $T B D$ is then substituted into the momentum corrector which calculates a new $T B D$ .

Since $T B D$ is ultimately recalculated by the momentum corrector — a momentum equation in explicit form — _is it necessary_ to solve the momentum predictor at all? The answer is “no”, unless it makes the solution more convergent which is more likely for high speed flows which are dominated more by momentum exchanges than mass conservation.

The transient algorithm should therefore include a simple switch that makes it possible to turn off the momentum predictor step. The rest of the algorithm remains the same, with the momentum matrix constructed to provide $T B D$ and $T B D$ .

With this switch included, the momentum predictor is then turned off in a lot of applications, _e.g._ highly viscous flows. In practice, it is often switched off for more complex flows, _e.g._ multiphase flows, which are beyond the scope of this book.



== The PIMPLE algorithm

The $T B D$ - $T B D$ coupling algorithms in Sec. 5.12 and Sec. 5.19 can be combined into an algorithm known as PIMPLE. PIMPLE merges the controls of PISO and SIMPLE (hence the merged acronym), in particular the iterative loops and under-relaxation.

All controls are optional; the standard transient algorithm is replicated by deactivating both the under-relaxation and the PIMPLE loop. By including the PIMPLE loop, equations are solved using variables updated within the time step. Accuracy is improved in particular due to the update of matrix coefficients from the contribution of $T B D$ to advection.

For transient simulations, temporal accuracy can be maintained at a higher $T B D$ ( $T B D$ ) using a second order time scheme (Sec. 3.18) combined with iterations of the PIMPLE loop. Similarly, the PIMPLE loop can update explicit source terms, _e.g._ in energy or momentum, to improve accuracy.

=== Pseudo-transient solution

PIMPLE can be configured to produce a steady flow solution quickly by a _pseudo-transient_ simulation. These simulations are not intended to capture realistic transient behaviour so can run at $T B D$ with some under-relaxation if necessary.

The simulations can be accelerated to a steady state using _local time stepping_ (LTS). LTS recognises that $T B D$ is limited by the maximum $T B D$ associated with the cell with small $T B D$ and/or high $T B D$ . It uses a field of $T B D$ corresponding to the $T B D$ limit _in each cell_ to accelerate the transient solution. While using a $T B D$ field makes the transient solution invalid, it is acceptable at steady state when $T B D$ .



== Solving for energy

Examples with energy conservation in this chapter assume $T B D$ is constant to produce the transport equation for $T B D$ in Eq. (2.65). When $T B D$ cannot be assumed to be constant, it is common to solve an equation for internal energy $T B D$ instead, _e.g._ Eq. (2.60) with the material derivative replaced using Eq. (2.14)

$ T B D $
This equation can be discretised to form a matrix equation $T B D$ , from which $T B D$ can be solved. The challenge is that the diffusion term $T B D$ is not expressed in terms of $T B D$ , so is discretised explicitly, which adversely affects convergence.
To improve convergence, an implicit term is introduced which is similar in form and scale to the problem term. For energy Eq. (5.39), we use $T B D$ , where the diffusivity $T B D$ is calculated from the $T B D$ function of $T B D$ .

The extra term is both “added and subtracted” in implicit and explicit form. This has the effect of increasing the contributions to matrix coefficients, while cancelling out a large part of the explicit contribution from $T B D$ , as illustrated below

$ T B D $
The overall solution procedure involves first updating $T B D$ from $T B D$ from thermodynamic relationships, _e.g._ Eq. (2.62). The equation for $T B D$ is then solved, and the subsequent solution for $T B D$ is converted back to $T B D$ .
If $T B D$ is expressed as a polynomial function of $T B D$ , Eq. (2.64), then $T B D$ is converted to $T B D$ using the analytical integral of Eq. (2.62)

$ T B D $

=== From temperature to energy

The conversion from $T B D$ to $T B D$ is more complex because $T B D$ cannot be made the subject of Eq. (5.40). Instead, it can be “inverted” using an iterative scheme such as the Newton-Raphson method,#footnote[Isaac Newton,_De analysi per aequationes numero terminorum infinitas_, 1669 and Joseph Raphson, _Analysis aequationum universalis_, 1690.] which for this problem is:

$ T B D $
$T B D$ is updated from the current $T B D$ using the evaluated $T B D$ and $T B D$ polynomials, Eq. (5.40) and Eq. (2.64), respectively. One iteration is often sufficient for the error $T B D$ to fall within a prescribed tolerance, but further iterations of Eq. (5.41) can be applied if necessary.
=== Boundary conditions

The boundary conditions for the energy equation are generally specified in terms of $T B D$ . But since they must be applied to the variable being solved, they must be reposed in terms of $T B D$ . A fixed value condition $T B D$ is converted to an equivalent condition $T B D$ , _e.g._ by Eq. (5.40). A fixed gradient condition $T B D$ for $T B D$ is converted to a fixed gradient $T B D$ for $T B D$ by

$ T B D $



== Summary of algorithms and solvers

- The finite volume method produces sparse matrices, Sec. 5.1.
- Symmetric matrices have the same coefficients across the diagonal and advection generally produces an asymmetric matrix.
- The matrix equations are solved using iterative methods that reduce a residual to a specified tolerance, Sec. 5.4.
- Diagonal dominance is a sufficient condition for convergence, see Eq. (5.9) and Sec. 5.5.
- Under-relaxation may be required to maintain diagonal dominance, Sec. 5.6, in particular for steady-state solutions.


=== Gauss-Seidel method

- Simple iterative method for solutions with cost $T B D$ , Sec. 5.2.
- Effective when the range of influence is short, _i.e._  $T B D$ in each cell is influenced by changes in close neighbouring cells.
- The symmetric variant and mesh renumbering may improve convergence, Sec. 5.8.


=== Conjugate gradient method

- A “descent” method which solves a minimisation problem for cases with a longer range of influence, Sec. 5.14.
- Solved iteratively, which requires preconditioning to make it efficient, Sec. 5.16.
- Symmetric matrix: PCG with DIC preconditioning.
- Asymmetric matrix: PBiCGStab with DILU preconditioning.


=== Geometric-algebraic multi-grid (GAMG)

- Well suited to equations with a long range of influence, in particular the pressure equation, Sec. 5.17.
- Simple algebraic agglomeration is corrected by a scaling based on minimisation for Laplacian-dominated equations, Sec. 5.18.


=== Systems of equations

- Systems of equations are solved iteratively with solutions from one equation substituted into the next, Sec. 5.9.
- Mass and momentum conservation are coupled by a pressure equation and momentum and flux correctors, Sec. 5.10.


=== Steady-state solution

- Steady flow solutions use the SIMPLE algorithm, Sec. 5.12.
- Equations are solved in an iterative sequence with under-relaxation required to ensure convergence, Sec. 5.13.


=== Transient solution

- Transient solutions use a PISO loop, Sec. 5.19.
- Correction for non-orthogonality in the pressure equation may require a distinct corrector loop, Sec. 5.20.
- In the PISO loop, typical configurations solve: 2 pressure equations, each including 1 non-orthogonal corrector; or 3 pressure equations without any non-orthogonal correctors.
- The PIMPLE algorithm provides additional controls for transient and pseudo-transient simulations, Sec. 5.21.

For CFD to be useful, it must calculate the flow, forces, heat transfer, _etc._ with sufficient accuracy. The accuracy of any calculation is always compromised by various approximations.

First, there are approximations from the numerical methods that are employed. Discretisation schemes, described in Chapter 3, affect accuracy which can be improved by reducing the size of cells and time step, at the cost of increased solution times.

The algorithms and solvers, described in Chapter 5, are a further source of numerical approximation. Linear solvers are generally iterative with tolerances that balance computational cost against accuracy of the solution. Algorithms may include loops with a fixed number of iterations that limits the overall solution convergence.

Boundary conditions, described in Chapter 4, are applied as part of the configuration of the CFD simulation. Invariably there is some level of approximation in the way the conditions represent the problem of interest. The approximation is physical, not numerical, and includes everyday examples such as:

- setting a fixed value condition with a uniform value, _e.g._ on velocity at an inlet;
- applying a zero gradient condition, _e.g._ to pressure at a wall boundary;
- imposing the direction of inflow, _e.g._ with the inlet-outlet-velocity condition, Sec. 4.15.


Beyond the numerical methods and boundary conditions, are the models that describe fluid dynamics. The governing equations and basic models in Chapter 2 are accurate for laminar flows with modest speeds and heat transfer. However, as the flow physics becomes more complex, the models become increasingly approximate.

It is clear from Chapter 6 that turbulence is highly complex, so turbulence modelling is often a major source of error in a CFD solution. It exemplifies the observation of statistician George Box that _all models are_ _wrong_.#footnote[George Box, Science and statistics, 1976.] He advises that we should seek “an economical description of natural phenomena” and “simple but evocative models”.

It is with that advice we present some of the Reynolds-averaged turbulence models commonly used in CFD. The models take the form of transport equations for turbulence fields, typically turbulent kinetic energy and dissipation rate.

The models are accompanied by advice on specification of the inlet and initial conditions for the fields. Particular attention is given to the modelling of the flow within turbulent boundary layers at wall boundaries.



== The k-epsilon turbulence model

The $T B D$ *model* (k-epsilon) was the first turbulence model to be widely adopted for a variety of flows in CFD.#footnote[William Jones and Brian Launder, _The prediction of laminarization with a two-equation model of_ _turbulence_, 1972.] It is one of a family of _two equation_ models which solve two transport equations, one usually for $T B D$ and one for another variable, often a dissipation rate. These models are the industry standard for CFD.

The $T B D$ model combines the transport equations for $T B D$ and $T B D$ , Eq. (6.28) and Eq. (6.32) respectively, which are reproduced below assuming $T B D$  = constant and replacing material derivatives in conservative form as in Eq. (2.26).

$ T B D $

$ T B D $

The standard model coefficients are

$ T B D $
The choice of coefficients in Eq. (7.3) originates from a large series of validation simulations of free shear flows.#footnote[Brian Launder, Alan Morse, Wolfgang Rodi and Brian Spalding, _Prediction of free shear flows — A_ _comparison of the performance of six turbulence models_, NASA Conference on Free Shear Flows, Langley, 1972.]
The model can be deployed in a steady flow solution, setting the local time derivatives to zero, _e.g._ $T B D$ . It can also form part of a transient solution to capture unsteady features such as vortex shedding.

The model uses a momentum equation in ensemble-averaged form, _e.g._ Eq. (6.26), using $T B D$ calculated from updated $T B D$ and $T B D$ according to Eq. (6.27).

The model equations are added to the end of the main loop in the steady and transient algorithms from Sec. 5.12 or Sec. 5.21, respectively. The figure below shows the additions to the transient algorithm, starting from the momentum corrector step.

Following the momentum and flux correctors, $T B D$ can be updated, _e.g._ for non-Newtonian models where $T B D$ . $T B D$ and $T B D$ are calculated and stored for source terms in both the $T B D$ and $T B D$ equations. Those equations are then solved, followed by an update to $T B D$ , and subsequently $T B D$ .



== Initialisation of the k-epsilon model

Initial values and boundary conditions must be specified for $T B D$ and $T B D$ to solve their respective transport equations. The ideal specification of boundary conditions for $T B D$ and $T B D$ follows those for $T B D$ described in Sec. 4.3.

Turbulence fields require: a fixed value condition at inlets; zero gradient or inlet-outlet at outlets; and, a more complex specification at solid walls, introduced in Sec. 7.7.

Inlet values $T B D$ and $T B D$ must therefore be specified. There may be industry standards, published recommendations or measured data to help select these values for the specific problem being simulated.

But more often than not, $T B D$ and $T B D$ must be estimated. Inlet and initial estimates of $T B D$ are usually based on a Prandtl mixing length $T B D$ from the expression

$ T B D $
This relation can be derived by considering a planar boundary layer. In the “log law” region (see Sec. 7.7) $T B D$ which can be combined with the mixing length Eq. (6.21) to give
$ T B D $
Substituting $T B D$ from Eq. (6.31) yields Eq. (7.4).
A value for $T B D$ must then be specified in order to calculate inlet and initial values of $T B D$ from Eq. (7.4). Procedures to estimate $T B D$ are described in Sec. 7.3.

Inlet and initial estimates for $T B D$ can be calculated by

$ T B D $
from a specified _turbulent intensity_ $T B D$ , the ratio of the root-mean-square (RMS) of turbulent fluctuations $T B D$ to the mean flow speed $T B D$ . The expression is derived from the definition $T B D$ .
A value for $T B D$ must then be specified in order to calculate the inlet and initial values of $T B D$ from Eq. (7.6). Procedures to estimate $T B D$ are also described in Sec. 7.3. The values of $T B D$ and $T B D$ at inlet boundaries influence the solution throughout the CFD simulation, so should be estimated as accurately as possible.

The _accuracy_ of the initial (internal) values is not so critical, since they do not influence the solution beyond a short period at the beginning of a simulation.

Initial values can, however, affect _stability_ during the early steps of a CFD simulation. The flow boundary conditions generally cause sudden impulses which can generate large forces, causing fluctuations in the solution of $T B D$ . Higher $T B D$ , based on initial $T B D$ and $T B D$ values, tends to cause larger fluctuations, which may make the solution of the $T B D$ -equation unstable.



== Inlet turbulence

Expressions are presented in Sec. 7.2 to estimate inlet and initial values of $T B D$ and $T B D$ . They include parameters $T B D$ and $T B D$ which must themselves be estimated sufficiently accurately to calculate $T B D$ and $T B D$ reliably.

The values of $T B D$ and $T B D$ at domain inlets depend on the flow conditions upstream of the inlet. The figure below shows typical ranges of intensity $T B D$ for different upstream flow conditions.

A medium intensity $T B D$ is most commonly specified in CFD problems, in particular for internal flows. For these flows, $T B D$ can be calculated from a power-law function of $T B D$ , fitted to measurements at the central axis in fully developed flow along a smooth-wall pipe, according to#footnote[Nils Basse, _Turbulence intensity and the friction factor for smooth- and rough-wall pipe flow_, 2017.]

$ T B D $
An estimate of $T B D$ at the centre axis of a pipe, see Sec. 6.12, can be used in conjunction with $T B D$ from of Eq. (7.7). For ducts and channels of non-circular cross-section, $T B D$ can be calculated by $T B D$ , where $T B D$ is the cross-sectional area and $T B D$ is the perimeter length. For a partially filled pipe or duct, $T B D$ corresponds to the wetted region where the fluid is in contact with the boundary.
For wall-bounded flows with a boundary layer of thickness $T B D$ , an estimate of $T B D$ is often used. This relation (see also Sec. 6.12) requires $T B D$ to be estimated, _e.g._ from the $T B D$ expression for a turbulent layer at the end of Sec. 6.4.

=== Verifying turbulent viscosity

Combining Eq. (7.4), Eq. (7.6) and Eq. (6.31) gives the following expression for $T B D$ in terms of length $T B D$ and velocity $T B D$ scales:

$ T B D $
Values of $T B D$ need to be realistic. Realistic values usually fall within the range of molecular viscosities $T B D$ for common fluids at standard temperature shown below.
The range is presented in terms of _kinematic_ viscosity $T B D$ which governs the _rate_ of momentum _diffusion_, _e.g._ the rate of growth of boundary layers. By contrast, _forces_ are governed by dynamic viscosity $T B D$ , which make liquids “feel” more viscous.



== Turbulent boundary layers

At solid walls, the tangential flow speed $T B D$ increases rapidly across a thin boundary layer, as discussed in Sec. 6.4. At high $T B D$ , the velocity profile has a universal character shown below.

The profile compares measured data#footnote[Joseph Kestin and Peter Richardson, _Heat transfer across turbulent, incompressible boundary layers_, 1963.] in terms of a dimensionless velocity $T B D$ and distance to the wall $T B D$ , given by

$ T B D $
Both parameters in Eq. (7.9) are based on a _friction velocity_ $T B D$ which is related to the wall shear stress $T B D$ by
$ T B D $
At the wall $T B D$ . Close to the wall, $T B D$ is suppressed, creating a region where flow is laminar $T B D$ , known as the _viscous sub-layer_ . The profile in this region is described by the relation
$ T B D $
Turbulence becomes significant through the _buffer layer_ which describes the region $T B D$ . Van Driest provides a model for the increase in mixing length through this region, by#footnote[Edward van Driest, _On turbulent flow near a wall_, 1956.]
$ T B D $
Finally, in the _inertial sub-layer_ for $T B D$ , flow is turbulent and the velocity profile is described by the _logarithmic law of the wall_, often abbreviated to simply the _log law_, according to#footnote[Alternatively written $T B D$ , where $T B D$ .]
$ T B D $
The equation includes Kármán’s constant $T B D$ and constant $T B D$ . For a smooth wall, $T B D$ #footnote[Hermann Schlichting and Klaus Gersten, _Boundary-layer theory_, 2017.] – 5.5 is commonly used. Both Eq. (7.11) and Eq. (7.13) can be derived assuming a constant shear stress across the profile, equating to $T B D$ at the wall. In the viscous sub-layer, the shear stress is laminar, so
$ T B D $
This equation integrates with a zero constant of integration to give $T B D$ , from which Eq. (7.11) is derived. In the inertial sub-layer, the shear stress is turbulent (laminar is negligible), so
$ T B D $
Assuming $T B D$ gives $T B D$ , which integrates to yield Eq. (7.13). In the inertial sub-layer, $T B D$ as described in Eq. (7.5), which combines with Eq. (7.15) and Eq. (6.31) to give
$ T B D $



== Wall functions

CFD simulations may be used to calculate the forces on solid bodies exerted by the fluid, _e.g._ in aerodynamics. The wall shear stress is then calculated according to $T B D$ . With turbulent boundary layers, the $T B D$ calculation requires cells with very small lengths normal to the wall to be accurate. The resulting mesh is inevitably large, which carries a high computational cost. The problem for CFD is how to calculate $T B D$ with sufficient accuracy, but at an affordable cost.

_Wall functions_ provide a solution to this problem by exploiting the universal character of the velocity distribution described in Sec. 7.4. They use the law of the wall Eq. (7.13) as a model to provide a reasonable prediction of $T B D$ from a relatively inaccurate $T B D$ calculation at the wall.

Wall functions use the _near-wall cell centre height_ $T B D$ , _i.e._ the distance to the wall from the centre P of each near-wall cell. Typically when using wall functions, $T B D$ should correspond to a $T B D$ within the typical range of applicability of the log law Eq. (7.13), _i.e._ $T B D$ .

With such a mesh, the calculated $T B D$ is then significantly lower than its true value. Wall function models compensate for the resulting error in the prediction of $T B D$ by increasing viscosity at the wall. The increase is applied to $T B D$ at the wall patch faces, which would otherwise be $T B D$ , corresponding to $T B D$ .

=== Standard wall function

The *standard* wall function for a “smooth” wall calculates $T B D$ for each patch face based on the near-wall $T B D$ . No adjustment is made to $T B D$ when $T B D$ corresponds to the viscous sub-layer. When $T B D$ corresponds to the inertial sub-layer, $T B D$ is calculated by

$ T B D $
The condition ( $T B D$ ) that determines whether $T B D$ lies within the inertial sub-layer corresponds to $T B D$ at the intersection of Eq. (7.11) and Eq. (7.13), calculated iteratively as
$ T B D $
$T B D$ is calculated in each near-wall cell according to:
$ T B D $
This expression is derived from Eq. (7.9) for $T B D$ and Eq. (7.16). The subscript $T B D$ denotes all properties are evaluated at P.
The wall function Eq. (7.17) is derived from the notion that $T B D$ is calculated numerically (assuming a stationary wall) by:

$ T B D $
At the same time, combining Eq. (7.9) and Eq. (7.10) gives
$ T B D $
Comparing Eq. (7.20) and Eq. (7.21) gives $T B D$ , which combines with Eq. (7.13) to yield Eq. (7.17).


== Alternative wall functions

The standard wall function described in Sec. 7.5 uses a function for $T B D$ that is discontinuous at $T B D$ , switching to $T B D$ for $T B D$ . A *continuous* wall function is available which evaluates $T B D$ as $T B D$ from a single equation describing the universal character of the velocity profile at high $T B D$ ,#footnote[Brian Spalding, _A single formula for the law of the wall_, 1961.]

$ T B D $
where $T B D$ . The equation combines Eq. (7.11) and Eq. (7.13), “disabling” Eq. (7.13) at low $T B D$ by subtracting low order terms from a polynomial expansion of $T B D$ .
The wall function is applied by solving Eq. (7.22) for $T B D$ from $T B D$ calculated from Eq. (7.9) using the near-wall cell centre height $T B D$ . The friction velocity is calculated by $T B D$ , where $T B D$ is the near-wall cell velocity. Finally, $T B D$ on the wall patch is calculated from a numerical interpretation of Eq. (7.20),

$ T B D $
where $T B D$ is the surface-normal velocity gradient. An iterative method is required to invert Eq. (7.22) and to accommodate other nonlinearities, _e.g._ $T B D$ is itself a function of $T B D$ .
=== Rough wall function

The standard wall function in Sec. 7.5 is applicable to smooth walls so does not account for surface roughness. Roughness is significant when the roughness “scale” $T B D$ #footnote[‘ $T B D$ ’ is often used to denote roughness, but we use ‘ $T B D$ ’ to avoid confusion with turbulent kinetic energy.] becomes larger than the thickness of the viscous sub-layer.

At higher surface roughness, turbulent eddies are generated near the wall at a scale $T B D$ , rather than $T B D$ . The viscous effects become negligible, causing the non-dimensionalised distance to become $T B D$ in the log law Eq. (7.13). To reflect this, Eq. (7.13) is modified to a form

$ T B D $
where $T B D$ a roughness function, dependent on $T B D$ . An intuitive model for $T B D$ is#footnote[Cyril Colebrook, _Turbulent flow in pipes, with particular reference to the transitional region between_ _smooth and rough wall laws_, 1939.]
$ T B D $
This *rough* wall function Eq. (7.24) reduces to Eq. (7.13) using the conventional $T B D$ definition of Eq. (7.9) at $T B D$ . As $T B D$ , it reduces to Eq. (7.13) using $T B D$ .
It is open to interpretation how to determine $T B D$ from roughness measurements of a surface. The parameter is sometimes split into $T B D$ , where $T B D$ is a measured _sand grain roughness height_ and $T B D$ is a coefficient that depends on the shape, consistency and packing of the roughness elements. Using that approach, values of $T B D$ often yield a good match between Eq. (7.24) and measured data.



== Turbulence near walls

The characteristic velocity distribution in turbulent boundary layers in Sec. 7.4, provided wall functions expressed as boundary conditions for $T B D$ in Sec. 7.5 and Sec. 7.6. Boundary conditions also need to be specified for turbulence fields at solid walls.

The turbulence generation $T B D$ influences the distribution of turbulence fields near a wall. At the wall, $T B D$ . In the inertial sub-layer, $T B D$ from Eq. (7.5) with: $T B D$ , obtained by combining Eq. (6.21) and Eq. (7.15); and, Eq. (6.24).

Since $T B D$ decreases with increasing $T B D$ in the inertial sub-layer, $T B D$ passes through a peak within the buffer layer (at $T B D$ ).

The peak in $T B D$ causes a similar peak in $T B D$ , shown in the following diagram. To the left of the peak, turbulent energy is transported back towards the wall by diffusion $T B D$ The profile of dissipation $T B D$ results from $T B D$ , obtained from Eq. (7.1). Very close to the wall ( $T B D$ ), diffusion is predominately molecular, such that it non-zero at the wall. The dissipation $T B D$ also has a non-zero value $T B D$ at the wall.

=== Wall functions and turbulence fields

When using a turbulence model, such as the $T B D$ model described in Sec. 7.1, boundary conditions must be specified for $T B D$ and $T B D$ at solid walls. The distribution of $T B D$ , non-dimensionalised as $T B D$ , close to the wall is shown below.

At the wall, $T B D$ but it rises quickly to a peak at $T B D$ before levelling off at $T B D$ as $T B D$ .

With wall functions, the height of the centre of each near-wall cell should correspond to $T B D$ within the range $T B D$ . Viewed at that scale, the $T B D$ profile appears flat. For $T B D$ , there is no such simple profile shape. These observations lead to the boundary conditions for $T B D$ and $T B D$ when using wall functions:

- *zero gradient* for $T B D$ ;
- *calculated near-wall cell value* for $T B D$ , according to:


$ T B D $
The calculation uses $T B D$ from the near-wall cell. The expression for $T B D$ is Eq. (7.4) with Eq. (6.24). The expression for $T B D$ uses the asymptotic condition $T B D$ for $T B D$ in Eq. (7.28).


== Resolving the viscous sub-layer

Wall functions were introduced in Sec. 7.5 to avoid the need for a large mesh with cells small enough to resolve the boundary layer into the viscous sub-layer. They provide a reasonable prediction of $T B D$ using the log law for the velocity distribution in the inertial sub-layer.

A CFD simulation may alternatively use a mesh with sufficiently thin cells to resolve the flow through the viscous sub-layer, _e.g._ with near-wall cell centre height corresponding to $T B D$ = 1, for a more accurate prediction of $T B D$ . If so, the turbulence model must then be able to function reliably in viscous flow regions.

Such models are usually described as “low Reynolds number”. The expression does not refer to the $T B D$ of the flow based on the characteristic scales of the problem, _e.g._ axial mean flow speed and diameter for a pipe. Instead it is a “turbulence” Reynolds number $T B D$ based on the scales of speed $T B D$ and size $T B D$ of turbulent eddies and can be defined as

$ T B D $
This definition is obtained from scale arguments introduced in Sec. 6.6, in which $T B D$ . Since $T B D$ represents fluctuations, $T B D$ . Combining these expressions into a Reynolds number yields Eq. (7.27).
=== Asymptotic consistency

Low- $T B D$ turbulence models pay attention to the behaviour of fluctuating velocities, _e.g._ $T B D$ , $T B D$ , in the limit that $T B D$ at the solid boundary.

They aim to capture the shape the profiles of $T B D$ and $T B D$ as they approach $T B D$ . Let $T B D$ and $T B D$ define the directions tangential and normal to the wall respectively. Profiles in the fluctuating velocities can be expressed by polynomials in $T B D$ , _i.e._

$ T B D $
where $T B D$ , _etc._ are functions of space and time. The no slip condition implies $T B D$ , so $T B D$ to the lowest order in $T B D$ . For $T B D$ , it is $T B D$ since, at the wall, $T B D$ and by continuity $T B D$ .
The turbulent properties are, to the lowest order in $T B D$ , as follows.

- From Sec. 6.11, $T B D$
- From Eq. (6.30), $T B D$


It follows that models achieve asymptotic consistency when

$ T B D $
in the limit that $T B D$ .


== Low-Re k-epsilon models

There are many low- $T B D$ turbulence models for CFD simulations where the cells near the solid walls are sufficiently thin to resolve the flow through the viscous sub-layer.

Among them are several low- $T B D$ $T B D$ models based on Eq. (7.1) and Eq. (7.2) with additional corrections $T B D$ , $T B D$ , $T B D$ and $T B D$ :

$ T B D $

$ T B D $
The calculation of $T B D$ also includes a correction $T B D$ :
$ T B D $
This model was first presented by *Jones and Launder* in their seminal $T B D$ publication.#footnote[William Jones and Brian Launder, _The prediction of laminarization with a two-equation model of_ _turbulence_, 1972.] They proposed functions for $T B D$ , $T B D$ , $T B D$ , $T B D$ and $T B D$ , as well as the coefficients $T B D$ , $T B D$ and $T B D$ .
*Launder and Sharma* subsequently presented#footnote[Brian Launder and B.I. Sharma, _Application of the energy-dissipation model of turbulence to the_ _calculation of flow near a spinning disc_, 1974.] the model with a modified $T B D$ function and the more established coefficients listed in Eq. (7.3).

The resulting model became known as the Launder-Sharma $T B D$ model.#footnote[although *Jones-Launder-Sharma* $T B D$ *model* would seem more equitable.] It is arguably the most popular low- $T B D$ model today.

The first notable modification to the standard $T B D$ model is $T B D$ (sometimes denoted “ $T B D$ ”) in Eq. (7.29). It is the dissipation rate at the wall ( $T B D$ ), see figure, Sec. 7.7, calculated by

$ T B D $
The term equates to $T B D$ in the boundary layer which is consistent with Eq. (7.28). The benefit of redefining the dissipation rate as $T B D$ is that the boundary conditions at a wall for the Launder-Sharma model are the same for $T B D$ and $T B D$ :
- *fixed value* $T B D$ ;
- *fixed value* $T B D$ .


The next significant modification is the $T B D$ function

$ T B D $
This modification recognises that $T B D$ so decreases $T B D$ through the buffer and viscous sub-layer to the wall, consistent with the decrease in $T B D$ according to the van Driest model Eq. (7.12).
The extra term $T B D$ in Eq. (7.30) is a follows, designed so that $T B D$ matches its recognised peak value within the buffer layer:

$ T B D $
Finally, $T B D$ and $T B D$ provide damping of the production and dissipation terms close to the wall in Eq. (7.30). The standard functions are $T B D$ (_i.e._ no damping), and
$ T B D $
which gives $T B D$ at the wall.


== Specific dissipation rate

The $T B D$ model is one of a family of _two-equation_ models for turbulence. With two equations, the models can represent each of the scales, $T B D$ and $T B D$ , which characterise $T B D$ . Most often, $T B D$ is used to represent $T B D$ .

The other variable must represent $T B D$ and so far we have used $T B D$ with SI units $T B D$ . The _specific dissipation rate_ $T B D$ , with SI units of $T B D$ , is a popular alternative for this variable in turbulence modelling.

While Kolmogorov first proposed a two-equation $T B D$ model,#footnote[Andrey Nikolaevich Kolmogorov, _Equations of turbulent motion in an incompressible fluid_, first published in Russian in _Izv. Akad. Nauk SSSR_ *6*, 1941.] the $T B D$ models used in CFD originate from Wilcox.#footnote[David Wilcox, _Reassessment of the scale-determining equation for advanced turbulence models_, 1988.] Here, “models” is plural since there are several versions of $T B D$ model with modifications and additions from its original form.

The *original* $T B D$ *model* is presented below (with some changes to the original variable names), assuming $T B D$  = constant for direct comparison with $T B D$ in Sec. 7.1.

$ T B D $

$ T B D $

The standard model coefficients are

$ T B D $
Comparing dissipation terms in Eq. (7.1) and Eq. (7.36) gives the relation $T B D$ . Substituting in Eq. (6.31) leads to a simple relation for turbulent viscosity, given by
$ T B D $
Inlet and initial estimates for $T B D$ can be calculated by
$ T B D $
using $T B D$ , in a manner similar to $T B D$ in Eq. (7.4).
With wall functions, the _boundary condition_ applied to $T B D$ set a *near-wall cell value* according to

$ T B D $
The expression for $T B D$ ( $T B D$ ) is a solution to the following equation for the viscous sub-layer where diffusion and dissipation terms dominate in Eq. (7.37):
$ T B D $
The equivalent dissipation terms for $T B D$ in Eq. (7.37) and $T B D$ in Eq. (7.2) are $T B D$ and $T B D$ respectively. The former is more stable in a numerical solution since it is insensitive to variations in $T B D$ .


== Enhancements to the k-omega model

A comparison of the $T B D$ and $T B D$ models shows the dissipation term is $T B D$ in Eq. (7.37) for $T B D$ and $T B D$ in Eq. (7.2) for $T B D$ . While the $T B D$ term behaves well, the $T B D$ term is troublesome at a wall as $T B D$ , so requires damping in the low- $T B D$ formulation in Sec. 7.9 to resolve the viscous sub-layer.

Its better dissipation term gives the $T B D$ model an advantage over $T B D$ in the near-wall region. The disadvantage the original $T B D$ model is its sensitivity to the freestream values of $T B D$ , which is not present in $T B D$ in the $T B D$ model. Neither model, in their original form, performs well under adverse pressure gradients.

However, since its initial publication, many enhancements have been made to the original $T B D$ model, in particular to address the problems mentioned above.

=== Cross-diffusion

The dependency on freestream values of $T B D$ is addressed by the inclusion of a _cross-diffusion_ term $T B D$ in the $T B D$ equation.#footnote[Florian Menter, _Zonal two equation_ $T B D$ _models for aerodynamic flows_, 1993.] The term is derived when the $T B D$ equation, _e.g._ Eq. (7.2), is expressed in terms of $T B D$ by substituting $T B D$ . Its form is due to the expansion of $T B D$ in the diffusion term by Eq. (2.74b).

The cross-diffusion term makes the $T B D$ equation more equivalent to $T B D$ , and thus independent of freestream values.

=== Stress limiter

The original $T B D$ and $T B D$ models are known to delay or suppress flow separation under adverse pressure gradients (described in Sec. 6.5). Under such conditions the ratio of the production to dissipation of $T B D$ can be significantly higher than unity. The calculated $T B D$ from Eq. (7.39) is excessively high, causing an over-prediction of shear stress $T B D$ .

The problem is alleviated by limiting the shear stress, based on the assumption it is proportional to $T B D$ in the boundary layer, _i.e._ $T B D$ where $T B D$ is a constant. A stress limiter is implemented through a modification to the calculation of $T B D$ :

$ T B D $
where $T B D$ is a 3D representation of $T B D$ .
=== Standard models

Different versions of the $T B D$ model are well catalogued in the _Turbulence modeling_ _resource_, NASA Langley Research Center, #link("https://turbmodels.larc.nasa.gov")[https://turbmodels.larc.nasa.gov].

Today, there are arguably two “standard” $T B D$ models, First, the $T B D$ *-* $T B D$ *SST* *model*,#footnote[Florian Menter, _Two-equation eddy-viscosity turbulence models for engineering applications_, 1994.] (SST = shear stress transport) which emerged as a popular choice in CFD over recent decades.

It combines the $T B D$ model near the wall with $T B D$ (expressed in terms of $T B D$ ) in the freestream, by applying _blending functions_ to model coefficients, the cross-diffusion term and the stress limiter.

Secondly, the $T B D$ *-* $T B D$ *2006 model*#footnote[David Wilcox, _Turbulence modeling for CFD, 3rd ed._, 2006.] applies the cross-diffusion term and stress limiting to the original $T B D$ model. The terms are applied using _switches_ so that the model maintains its simplicity, without the need for blending functions.



== Heat transfer in turbulent flow

The initial focus of turbulence modelling is to capture the effect of mixing on momentum diffusion since it influences the overall flow solution. But other properties are also transported by the turbulent eddying motions, in particular heat.

The effects of turbulence on heat transfer can be described using the following equation for internal energy $T B D$ , obtained by substituting the material derivative in Eq. (2.57) and ignoring $T B D$ :

$ T B D $
In turbulent flow, internal energy can be decomposed into averaged and fluctuating components $T B D$ , see Eq. (6.11). The ensemble average of the terms in $T B D$ introduces a heat flux
$ T B D $
This additional heat flux in the energy equation is equivalent to the Reynolds stress $T B D$ , Eq. (6.14), in momentum. Boussinesq modelled $T B D$ by Eq. (6.20), using the concept of an eddy viscosity $T B D$ corresponding to turbulent mixing, analogous to viscosity due to molecular motion according to Newton’s law Eq. (2.41).
Similarly, $T B D$ can be modelled using a turbulent thermal conductivity $T B D$ due to turbulent mixing, by analogy with Fourier’s law Eq. (2.54) for conduction due to molecular interaction

$ T B D $
The total heat flux $T B D$ in Eq. (7.44) is then expressed in terms of the combined turbulent mixing and molecular interaction, using an _effective thermal conductivity_ $T B D$ , as follows:
$ T B D $

=== Modelling turbulent heat transfer

Turbulent heat transfer can be incorporated into turbulence models based on eddy-viscosity and Reynolds-averaging, with additional thermal wall functions.

First, the calculation of $T B D$ by Eq. (7.47) requires $T B D$ from the turbulence model. A common approach to calculate $T B D$ is from $T B D$ based on an estimate of _turbulent Prandtl_ _number_

$ T B D $
$T B D$ provides a good estimate for many fluids, with $T B D$ often chosen as a default value for CFD calculations. For some more unusual fluids, _e.g._ liquid metals, $T B D$ .
=== Wall heat flux

The calculation of heat transfer through boundary walls is an important aspect of a many CFD simulations. Near walls, the distribution of $T B D$ tends to mimic $T B D$ .

Consequently, the challenges of calculating wall heat flux $T B D$ are similar to wall shear stress $T B D$ . Cells close to the wall must be very thin to resolve the viscous sub-layer in $T B D$ (when $T B D$ ).

Otherwise, wall functions can be used to adjust $T B D$ to compensate for the under-prediction of $T B D$ as described in Sec. 7.14.



== Thermal boundary layers

In a turbulent boundary layer, the distribution of temperature is similar to velocity, with the viscous and intertial sub-layers, separated by a buffer layer, as discussed in Sec. 7.4.

By analogy with $T B D$ , Eq. (7.10), we define friction temperature as

$ T B D $
The wall layer is then described by a dimensionless temperature
$ T B D $
where $T B D$ is the fluid temperature at the wall. Ignoring heat generation by viscous stresses, the profile in the viscous sub-layer is described by the relation
$ T B D $
The profile in the inertial sub-layer is commonly described by a log law for $T B D$
$ T B D $
The derivation of Eq. (7.51) and Eq. (7.52) assumes a constant heat flux across the profile, equating to $T B D$ at the wall. In the viscous sub-layer, the heat flux is laminar so $T B D$ and
$ T B D $
This equation integrates between $T B D$ at a distance $T B D$ from the wall to $T B D$ at the wall, to yield Eq. (7.51). In the inertial layer, the heat flux is assumed turbulent $T B D$ and
$ T B D $
Combining Eq. (6.24), Eq. (7.9) and Eq. (7.15) yields the ratio $T B D$ . Substituting in Eq. (7.54) and integrating then leads to Eq. (7.52) where $T B D$ is the constant of integration.
The constant $T B D$ is generally considered to be a function of $T B D$ . A reasonable approximation for this function is#footnote[Hermann Schlichting and Klaus Gersten, _Boundary-layer theory_, 2017.]

$ T B D $
Another function, commonly used in thermal wall functions, is $T B D$ , where $T B D$ is the function of $T B D$ :#footnote[Chandra Jayatilleke, _The influence of Prandtl number and surface roughness on the resistance of the_ _laminar sub-layer to momentum and heat transfer_, 1966.]
$ T B D $
The expression for $T B D$ uses the coefficient $T B D$ from Eq. (7.11). These constants of integration are sometimes subsumed within the log function as a coefficient “ $T B D$ ” in the log law expressions, as the footnote on page 483 explains.


== Thermal wall functions

Wall functions were introduced in Sec. 7.5 in order to improve the calculation of wall shear stress $T B D$ when cells are too large near a wall to resolve $T B D$ accurately. The same problem exists with heat flux $T B D$ and an under-predicted $T B D$ . As before, the universal character of the boundary layer can be exploited, this time to improve the calculation of $T B D$ .

The temperature distribution is characterised by Eq. (7.51) for the viscous sub-layer, and the log law Eq. (7.52) for the inertial sub-layer. The transition $T B D$ for $T B D$ occurs at the intersection of the two equations, _i.e._ when

$ T B D $
While the transition between these regimes is fixed at $T B D$ by Eq. (7.18) for the $T B D$ profile, the (iterative) solution of Eq. (7.57) is dependent on $T B D$ and $T B D$ .
Using $T B D$ and Eq. (7.55) for $T B D$ , $T B D$ for air at $T B D$ with $T B D$ . For water under the same conditions, $T B D$ and the corresponding $T B D$ .

A wall function can be derived which adjusts the turbulent conductivity $T B D$ , in a similar manner to $T B D$ in the standard wall function in Sec. 7.5. The model calculates $T B D$ for each patch face based on the near-wall cell $T B D$ .

No adjustment is made to $T B D$ when $T B D$ corresponds to the viscous sub-layer. When $T B D$ corresponds to the inertial sub-layer, $T B D$ is calculated as

$ T B D $
where $T B D$ denotes the laminar thermal conductivity, _i.e._ $T B D$ from Eq. (2.54), to distinguish it from Kármán’s constant $T B D$ . Eq. (7.58) uses $T B D$ from Eq. (7.19), as in the standard $T B D$ wall function.
The wall function is derived based on adjusting $T B D$ to improve the numerical calculation of $T B D$ by

$ T B D $
where $T B D$ represents a value close to the wall, _e.g._ in a near-wall cell. By comparison, Eq. (7.49) and Eq. (7.50) combine to give
$ T B D $
The heat flux $T B D$ is consistent between Eq. (7.59) and Eq. (7.60) when
$ T B D $
When the log law Eq. (7.52) is then substituted in Eq. (7.61), it provides the thermal wall function model which adjusts $T B D$ according to Eq. (7.58).


== Summary of turbulence modelling

- The $T B D$ model solves transport equations for turbulent kinetic energy $T B D$ and dissipation rate $T B D$ , Sec. 7.1.
- It is the original of a family of two-equation models, which ultimately provide $T B D$ to calculate the turbulent stresses.
- Initial and inlet values for $T B D$ and $T B D$ must be specified, which can be calculated from turbulent intensity $T B D$ and mixing length scale $T B D$ , respectively, Sec. 7.2.
- $T B D$ and $T B D$ are often estimated using functions that fit experimental data for fully developed turbulent flow, Sec. 7.3.
- The $T B D$ models replace $T B D$ by specific dissipation rate $T B D$ , with equivalent expressions for initial and inlet values, Sec. 7.10.
- The two “standard” $T B D$ models today are the $T B D$ SST model and $T B D$ 2006 model, Sec. 7.11.
- Models can provide a turbulent conductivity $T B D$ to calculate the turbulent heat flux for thermal problems, Sec. 7.12.


=== Turbulent boundary layers

- Turbulent boundary layers include a thin viscous sub-layer adjacent to the boundary with a linear velocity profile and, further from the boundary, the inertial sub-layer with a log law profile, Sec. 7.4.
- Profiles in temperature are similar to those for velocity, with equivalent linear and log law relationships, Sec. 7.13.
- Very thin cells are generally needed to resolve the viscous sub-layer to calculate the velocity gradient $T B D$ at the wall accurately, Sec. 7.5.
- Such thin cells within the boundary layer region can increase the mesh to a size which is prohibitively costly to run.


=== Wall functions

- Wall functions permit much larger cells near the wall, by exploiting the universal character of the velocity distribution.
- The functions increase $T B D$ at the wall to compensate for the under-prediction of $T B D$ with larger cells, to improve the prediction of the wall shear stress, Sec. 7.5.
- Thermal wall functions similarly increase $T B D$ at the wall to compensate for the under-prediction of $T B D$ , in order to improve the prediction of the wall heat flux, Sec. 7.14.
- Standard wall functions make no adjustment to $T B D$ when the near wall cell centre falls below the transition within a buffer layer, Sec. 7.5.
- Other models include a continuous function of $T B D$ through the viscous sub-layer to the wall and adjustments for surface roughness, Sec. 7.6.
- Boundary conditions for turbulence fields with wall functions are based on observed profiles of those fields, Sec. 7.7.


=== Models with resolved boundary layers

- Turbulence models must predict the universal character of boundary layers when the viscous sub-layer is resolved with sufficiently thin cells, Sec. 7.8.
- Models like the Launder-Sharma $T B D$ include source terms and damping functions to improve the predictions and to simplify boundary conditions, Sec. 7.9.

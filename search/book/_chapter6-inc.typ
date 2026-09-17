The previous chapters provide some general principles for CFD, from the governing equations and common models, to numerical methods, boundary conditions, algorithms and solvers.

Since the book is about general principles, it does not generally cover complex physical modelling. The exception is _turbulence_ which need to be addressed since it is so common in fluid flow and therefore a critical part of most CFD simulations.

Turbulence is illustrated by a plume of smoke below. The smoke rises vertically from the base, initially tracing a relatively straight, narrow path, characteristic of _laminar flow_. At some point, there is an abrupt transition to _turbulent flow_ in which fluid particles follow curved paths which cross one another in a disorderly, irregular manner, causing $T B D$ and $T B D$ to fluctuate randomly in time.

Turbulent flow evidently involves significant _mixing_ of the fluid which is an important consideration in many CFD problems. For example, it directly affects the rate of pollutant dispersion and the rate of reaction of a chemical process. It causes heat to diffuse more rapidly. Similarly, mixing affects momentum diffusion and, consequently, the aerodynamic forces within a fluid and surrounding vehicles, wind turbines, buildings, _etc._

Turbulence is a vast subject with an abundance of textbooks for various levels of expertise. As an introduction to turbulence, this chapter aims to provide:

- a short description of the nature of turbulence;
- an understanding of its importance in analysis and particular challenges it presents in CFD;
- the basic concepts of turbulence modelling to lay the groundwork for the models presented in Chapter 7.




== Reynolds experiment

Reynolds distinguished the laminar and turbulent flow regimes in his influential experiments in the early 1880s.#footnote[Osborne Reynolds, _An experimental investigation of the circumstances which determine whether_ _the motion of water shall be direct or sinuous, and of the law of resistance in parallel channels_, 1883.] At that time, Poisueille’s law predicted correctly that the resistance to flow in a pipe $T B D$ for laminar flow at speed $T B D$ , which typically occurs at small $T B D$ and/or diameter $T B D$ . Reynolds wanted to understand why, at higher $T B D$ and/or $T B D$ , resistance $T B D$ (approximately).

His experiment used a large, glass-walled tank filled with water. A glass tube with a flared intake passed through the tank and out through one wall. Water from the tank flowed along the tube at speed $T B D$ , controlled by an outlet valve. A jet of liquid dye was injected at the inlet to the tube to visualise the flow.

The observed behaviour is shown in Reynolds’s original drawings on the following page. At sufficiently low speed, the streak of dye followed a straight line along the tube, indicating laminar flow. The speed was increased in small steps until at some distance from the tube intake (typically, 30 $T B D$ ) the dye would mix with the water, rapidly filling the tube. This marked the _transition to turbulent flow_.

Reynolds made the following important observations.

- Distinct _eddies_ were visible during turbulent flow, see above.
- The transition speed was _very sensitive to disturbances_ in the water entering the tube, even due to the temperature variations in the water.
- Before full transition was reached, intermittent “spots” of turbulence appeared and disappeared.


=== Reynolds number

Reynolds argued that the transition is controlled by the Reynolds number, $T B D$ from Eq. (2.68), using the characteristic length $T B D$ in a pipe. His deduction was based on the idea of scale similarity, introduced in Sec. 2.21.

He observed that with minimal disturbance in the tank, the transition to turbulent flow occurred at $T B D$  13000. When disturbances were evidently present, the transition occurred at $T B D$  2000. An updated view for pipe flow is that for $T B D$  2000, all disturbances will decay preventing the onset of turbulent flow. For $T B D$  2000, transition depends on initial disturbances and the roughness of the pipe wall.



== A picture of turbulence

In turbulent flow, the fluid follows irregular curved paths, known as eddies. Eddies are generally illustrated by curved arrows, which does little to reveal the true nature of turbulence.

Instead, we can build a picture of turbulence from: photographs and videos from experiments, examples from engineering and everyday experiences; and, from images and animations from detailed computer simulations.

Generally, we see a mass of intertwined, flow structures of many different sizes. The structures move and rotate continuously. Through interactions with one another and the main flow field, they change shape and size rapidly.

The structures indicate regions of _vorticity_ $T B D$ in the fluid, as defined in Sec. 2.11. It is vorticity, not velocity, which is often used to describe turbulence. For example, Davidson#footnote[Peter Davidson, _Turbulence: an introduction for scientists and engineers_, 2004.] describes _turbulence_ as “a spatially complex distribution of vorticity which advects itself in a chaotic manner”, and a _turbulent eddy_ as “a blob, tube or sheet of vorticity and the associated flow.”

=== Vorticity generation

Vorticity is primarily generated at solid boundaries and can be categorised as “fast” or “slow”.#footnote[Bruce Morton, _The generation and decay of vorticity_, 1984.] Fast generation is localised to regions of relative acceleration between the solid and fluid, _e.g._ at the leading edge in flow over a semi-infinite flat plate or at the lip at the exhaust of the air cannon in Sec. 2.11.

Slow generation occurs over longer sections of boundary surface, _e.g._ the wall of a pipe with internal fluid flow. It is caused by the tangential pressure gradient at the surface.

Since vorticity is mostly generated at solid boundaries, we can try to imagine the how vortex structures might evolve. Consider a layer of a viscoelastic fluid on a flat surface — something like a sticky adhesive or putty. When shear is applied, the material tends to rotate rather than slide horizonatally.

Sections of material can break off and form tubes which roll across the surface. They tend to extend along their length and bend under uneven shear.

When we apply shear by hand, we can feel that the rotation is initiated within the uneven points, _e.g._ around the finger joints. In a similar way, the onset of turbulence is sensitive to surface roughness with the transition from laminar to turbulent flow occurring at higher $T B D$ for smoother surfaces.



== Vorticity transport

For an incompressible fluid (with other simplifying assumptions), vorticity obeys the transport equation , Eq. (6.2). This is a typical advection-diffusion equation, similar to Eq. (2.65) for heat, which is expressed in terms of the local time derivative and advection in conservative form by

$ T B D $
Notably, Eq. (6.1) does not include a term in $T B D$ . This is in contrast with conservation of (linear) momentum which can redistribute a perturbation in $T B D$ instantaneously across all the domain through $T B D$ , as discussed in Sec. 2.22.
Instead, like heat, vorticity evolves locally _only_, with a range of influence limited by advective and diffusive transport, as discussed on page 126.

Advection of vorticity is clearly illustrated by the smoke ring shown in Sec. 2.11. Diffusion occurs by viscous torques transferring angular momentum between fluid elements.

The source of vorticity, $T B D$ is due to vortices changing shape under the influence of a velocity gradient $T B D$ . If a vortex is stretched, _e.g._ under shear as shown above, its radius decreases, so angular velocity, and thus $T B D$ , increases. Similarly, $T B D$ decreases if the vortex is compressed.

=== The vorticity transport equation

For incompressible flow with constant $T B D$ (and zero, or constant, body force), vorticity obeys the following transport equation:

$ T B D $
The equation is derived from momentum conservation of a homogeneous, incompressible, Newtonian fluid with $T B D$  constant and ignoring body forces, Eq. (2.49). Noting that $T B D$ due to Eq. (2.46) leads to
$ T B D $
The vorticity equation is _derived from the curl (_ $T B D$ _) of Eq.__(__6.3__)_ combined with the vorticity definition Eq. (2.37). The first term is $T B D$ .
Replacing $T B D$ and $T B D$ by $T B D$ in Eq. (2.72d) and applying Eq. (2.37) gives

$ T B D $
The curl of the second term in Eq. (6.3) is expressed by the curl of Eq. (6.4) with Eq. (2.75a), Eq. (2.73a), Eq. (2.75e) and Eq. (2.46), leading to
$ T B D $
The curl of the third term $T B D$ by Eq. (2.74f). The curl of the fourth term $T B D$ by Eq. (2.75a).
Combining all the terms and applying the material derivative Eq. (2.14), $T B D$ to Eq. (6.5), leads to Eq. (6.2).



== Boundary layers

Boundary layers#footnote[Ludwig Prandtl, _Über Fl__üssigkeitsbewegung bei sehr kleiner Reibung_, 1904.] are regions of fluid formed along solid boundaries in which the velocity $T B D$ varies: from zero at the boundary (no-slip condition, Sec. 4.4); to a value largely unaffected by the proximity of the boundary, determined by the flow conditions.

The figure above shows a boundary layer for flow in the $T B D$ -direction at speed $T B D$ , along a flat solid boundary oriented in the $T B D$ -normal direction. At the boundary surface, the vorticity $T B D$ is significant.

Vorticity can be shown over a planar section of the boundary layer of width $T B D$ and height $T B D$ (above, right). Applying Stokes’s theorem, Eq. (2.39), the integral = $T B D$ along the upper line (with zero along the wall and the verticals sides). The average vorticity over plane area $T B D$ is therefore $T B D$ .

Boundary layers are the main source of vorticity for turbulence. Turbulence occurs when instabilities, _e.g._ induced by roughness of the boundary surface, cause the vorticity to become chaotic, sustained by a sufficiently high $T B D$ .

The growth of boundary layers is related to vorticity transport. For flow over a flat plate, vorticity generated at the leading edge is advected by the flow, while diffusing away from the plate.   \   \
Vorticity propagates by diffusion by a distance $T B D$ in time $T B D$ , see Sec. 2.22. In that time, it is advected a distance $T B D$ , where $T B D$ is the freestream flow speed. Comparing the distances over the same $T B D$ , the boundary layer thickness is
$ T B D $
The relation $T B D$ is suitable for _laminar_ boundary layers, with coefficient $T B D$ depending on the definition of $T B D$ . Data and analysis, including _e.g._ the Blasius solution#footnote[Paul Richard Heinrich Blasius, _Grenzschichten in Fl__üssigkeiten mit kleiner Reibung_, 1908.], indicate $T B D$ in the case of the “99% thickness”, _i.e._ the distance from the wall where velocity reaches 99% of its asymptotic value.
In _turbulent_ boundary layers, the diffusion front advances more rapidly due to mixing, see Sec. 6.11. As a result, $T B D$ is relatively insensitive to $T B D$ , _e.g._ the analytical solution $T B D$ , based on a one-seventh ( $T B D$ ) power law for the velocity profile.



== Boundary layer separation

Boundary layers provide the main source of vorticity for turbulence as discussed in Sec. 6.4. A boundary layer breaks away from the boundary when it reaches the end of the surface or by _separation_ before reaching the end. Vorticity and turbulence are thereby swept into regions of fluid away from solid boundaries.

Flow around a cylinder illustrates boundary layer separation, vorticity and turbulence. The fluid flows at uniform velocity $T B D$ upstream of the cylinder. It decelerates to _stagnation_, $T B D$ , at point A on the surface (in the $T B D$ -normal direction).

High pressure at A introduces a _favourable_ pressure gradient $T B D$ which increases the flow speed $T B D$ around the cylinder towards B, developing a boundary layer in the process.

The flow reaches a peak speed at B, then decelerates over the downstream side of the cylinder. The _adverse_ pressure gradient $T B D$ causes $T B D$ to decrease. At some point C, the velocity gradient can reach $T B D$ . Beyond C, the boundary layer can _separate_ such that $T B D$ along its profile, see point D.

Boundary layer separation in a cylinder depends on Reynolds number $T B D$ , Eq. (2.68), using $T B D$ and $T B D$ . For $T B D$ , there is no separation, with the flow exhibiting a pattern downstream that mirrors the upstream flow.#footnote[Sadatoshi Taneda, _Experimental investigation of the wakes behind cylinders and plates at low Reynolds_ _numbers_, 1956.]

For $T B D$ , the boundary layer separates with its vorticity sustaining a pair of vortices attached to the rear of the cylinder.

At $T B D$ , vorticity is released downstream as vortices break off from the cylinder in a periodic manner known as the _K__árm__án vortex street_, shown above.

At $T B D$ , the vorticity starts to become become chaotic, with turbulence beginning to appear in the vortices. At $T B D$ , the entire wake region becomes turbulent.

The frequency of vortex shedding is characterised by another dimensionless number from Eq. (2.68), the Strouhal number $T B D$ . For $T B D$ , experiments show#footnote[Vincenc Strouhal, _Über eine besondere Art der Tonerregung_, 1878 and Anatol Roshko, _On the_ _development of turbulent wakes from vortex streets_, NACA Report 1191, 1954.] $T B D$ , where $T B D$ is the period at which the vortex pattern repeats.



== Scales of turbulence

In Sec. 6.2, turbulence was pictured as a mass of intertwined eddies which are rotating and distorting. The rotation, _i.e._ vorticity, mostly originates at solid boundaries and propagates locally (Sec. 6.3) through boundary layers (Sec. 6.4), which can separate and shed into the flow field (Sec. 6.5).

Flow disturbances trigger instabilities which cause vortices to stretch, compress and/or ‘break up’. At a sufficiently high $T B D$ , coherent flow structures rapidly disintegrate into a mass of turbulent eddies. This process is described in Richardson’s poem:#footnote[Lewis Fry Richardson, _Weather prediction by numerical process_, 1922.]

#quote(block: true)[ Big whorls have little whorls that feed off their velocity, and little whorls have lesser whorls and so on to viscosity — in the molecular sense. ]

In other words, large eddies created by instabilities become progressively smaller until they reach a size where the dissipation of their kinetic energy due to (molecular) viscosity is significant. The loss of kinetic energy causes these eddies to “die out” before they can become any smaller.

=== Kolmogorov microscales

The _Kolmogorov microscales_#footnote[Andrey Nikolaevich Kolmogorov, _The local structure of turbulence in incompressible viscous_ _fluid for very large Reynolds numbers_, first published in Russian in _Dokl. Akad. Nauk SSSR_ *30*, 1941.] describe the smallest scales that can exist in turbulent flow. They can be derived from heuristic arguments as follows. Energy dissipates as heat at a rate of $T B D$ per unit mass, obtained from Eq. (2.58) and switching from ‘per unit volume’ by replacing $T B D$ by $T B D$ .

Using scale similarity arguments (Sec. 2.21) the average rate of dissipation of energy $T B D$ can then be estimated for the smallest scales, characterised by length $T B D$ and speed $T B D$ , to be

$ T B D $
where “ $T B D$ ” means ”of the order of magnitude”.
Dissipation occurs when the Reynolds number, based on the Kolmogorov scales, is of the order unity, _i.e._ $T B D$ . Combined with Eq. (6.7) this yields the Kolmogorov microscales below.

$ T B D $

At the largest scales, eddies generated from the mean flow can be characterised by a length $T B D$ and flow speed $T B D$ (of the scale of turbulent fluctuations). The kinetic energy per unit mass of these eddies $T B D$ .

Experiments show the lifetime of an eddy is of the order of the time for one revolution, the _turn-over_ time $T B D$ . Therefore the rate of transfer of kinetic energy from larger to smaller eddies is $T B D$ .

This rate of transfer of energy between scales must match the rate of dissipation, Eq. (6.7). Equating the two energy rates yields the correspondence between the Kolmogorov scales and the largest scales in terms of a large scale $T B D$ :

$ T B D $



== Energy cascade

The process of large eddies becoming smaller is introduced in Sec. 6.6. The process involves a transfer of kinetic energy from larger to smaller eddies, known as a the _energy_ _cascade_.

The variation in kinetic energy with eddy size is often illustrated by the energy spectrum graph above. The horizontal axis is wave number $T B D$ ,#footnote[The symbol $T B D$ , is commonly used in science to denote wave number; it should not be confused with its later use to denote turbulent kinetic energy.] which represents the number of eddies per unit length, _i.e._ $T B D$ where $T B D$  = eddy size.

The cascade starts with the largest eddies of size $T B D$ with the highest level of turbulent kinetic energy per unit mass (TKE). In the range, $T B D$ , TKE = $T B D$ .

Large eddies are usually _anisotropic_ due to the way turbulence is generated. For example, flow past a cylinder causes shedding of vortices whose shape are similar to the cylinder, _i.e._ longer in the direction of the cylinder axis.

Once moving through open space, eddies quickly stretch, bend, rotate and break, so that quite soon they become “blobs” of vorticity. The turn-over time, initially associated with shedding, decreases with eddy-size — by a factor of $T B D$ at the smallest scales, according to Eq. (6.9). Any dominant frequencies, _e.g._ those due to shedding, are thereby quickly lost.

Kolmogorov’s hypothesis was that the turbulence then becomes _isotropic_ and that $T B D$ is then only a function of $T B D$ and $T B D$ when $T B D$ is large. From this hypothesis emerges his _five-thirds law_ for the _inertial subrange_ $T B D$

$ T B D $
where $T B D$ for $T B D$ . The exponents of $T B D$ and $T B D$ are chosen to match the dimensions of $T B D$ which are $T B D$ , and noting $T B D$ .
Experimental data generally supports Eq. (6.10) — remarkably, given its simplicity.



== The cost of simulating turbulence

The Kolmogorov microscales are the smallest scales of length, velocity and time in turbulence. The fact that the scales are so small has important consequences for CFD.

In automotive external aerodynamics, boundary layers develop from the front of the vehicle which separate towards the rear. The flow recirculates, forming a wake region with high levels of turbulence.

The Kolmogorov length scale $T B D$ , with a time scale of $T B D$ , based on an estimated turbulent energy dissipation rate $T B D$ and $T B D$ for air.

To resolve the smallest eddies with CFD, a mesh with typical cell length $T B D$  mm would then be required. Since we would need to capture a wake region of volume $T B D$ , the corresponding mesh would consist of $T B D$ cells.

In automotive external aerodynamics today (2022), some of the largest simulations run with meshes of $T B D$ cells on $T B D$ processor cores. These mesh sizes are clearly several orders of magnitude smaller than those required for _direct numerical simulation_ (DNS) of turbulence at the smallest scales.

An important consequence of using a mesh that does not resolve the smallest scales is that the main mechanism of dissipating turbulent energy is not captured. This affects the rate of energy transfer within the energy cascade.

If the energy dissipation is neither captured nor introduced by diffusive numerical schemes, excess energy inevitably “backs up”, causing over-prediction of the turbulence in the resolved eddies. Similarly, turbulence is under-predicted with diffusive schemes that introduce excessive numerical dissipation.

_Large-eddy simulation_ (LES) instead uses modelling for the unresolved scales while capturing the large eddies using accurate numerical schemes. Effective models dissipate the energy in the unresolved scales at a realistic rate, _e.g._ by the five-thirds law.

Turbulent fluctuations are predicted within the resolved range of wave numbers, _i.e._ corresponding to lower frequencies. LES is therefore useful for problems involving excitation at those frequencies, _e.g._ flow-induced vibration, jet breakup and noise.

The simulations are still relatively costly, however, since they are inherently transient, requiring strong convergence within each time step for accuracy. Mesh sizes are still relatively large and long simulation times are needed to generate reliable time-averaged properties.

We therefore need a faster method to calculate turbulence, suitable for flows which are steady (or reasonably steady).



== Reynolds-averaged simulation

The computational cost of DNS and, to a lesser extent, LES is too great for most practical CFD, as discussed in Sec. 6.8. Instead, a _Reynolds-averaged simulation_ (RAS) provides a much more affordable method to calculate turbulence.

It solves equations for “averaged” field variables to avoid resolving small fluctuations. Rather than consider an average over a time interval, we imagine the same flow repeated several times under nominally the same initial conditions (2 examples below).

Solutions vary due of differences in initial conditions and the chaotic nature of turbulence. The _ensemble_ average calculates the mean solution $T B D$ for multiple _realisations_ of the same flow.

Each field, _e.g._ $T B D$ , is decomposed into the averaged field $T B D$ and field of random fluctuations $T B D$ , according to

$ T B D $
We can apply this decomposition to fields in the momentum conservation Eq. (2.27) with stress $T B D$ split into $T B D$ and $T B D$ as follows:
$ T B D $
Let us assume constant $T B D$ and that the body force $T B D$ is not subject to turbulent fluctuations. Splitting the remaining fields, _i.e._ $T B D$ , $T B D$ and $T B D$ , into instantaneous and fluctuating components according to Eq. (6.11), and taking the ensemble average of Eq. (6.12), yields the following:
$ T B D $
The averaged Eq. (6.13) is derived using the following relations for general fields $T B D$ and $T B D$ : $T B D$ ; $T B D$ ; $T B D$ ; $T B D$ . Relations for averaged derivatives are: $T B D$ ; $T B D$ ; $T B D$ .
=== Reynolds stress

The terms for mean quantities in Eq. (6.13) are the same as in Eq. (6.12). The difference is that Eq. (6.13) includes the additional $T B D$ term containing fluctuations $T B D$ .

The argument of this divergence derivative is a tensor known as the _Reynolds_ _stress_#footnote[Osborne Reynolds, _On the dynamical theory of incompressible viscous fluids and the determination of_ _the criterion_, 1895.]

$ T B D $
Substituting Eq. (6.14) in Eq. (6.13) eliminates fluctuation terms. The remaining equation is in terms of averaged properties only, so we can dispose of the average notation ( $T B D$ ) to give
$ T B D $
This ensemble-averaged equation is the same form as Eq. (6.13) but with the addition of $T B D$ . Solving this Reynolds-averaged equation is the key to low cost CFD with turbulence — but it requires a _model_ for the additional unknown $T B D$ .


== The nature of viscosity

Before discussing turbulence further, it is useful to examine the origins of viscosity. Viscosity is introduced in Sec. 2.12 of this book as part of the Newtonian constitutive model. The model was originally phenomenological, but was later derived directly from kinetic theory which describes a gas as a number of submicroscopic _particles_, _e.g._ atoms or molecules, in random motion.

The kinetic view of viscosity imagines a fluid in two dimensions, $T B D$ and $T B D$ , subjected to a shear force in the $T B D$ -direction. Although the mean flow is in the $T B D$ -direction, particles move in the $T B D$ -direction due to random fluctuations with a _mean speed_ $T B D$ .

Consider a plane at $T B D$ . A particle will pass through the plane if its path towards it is not interrupted by a collision which sends it moving away from the plane. Particles passing through the plane arrive from an average distance $T B D$ , where $T B D$ is some factor of the average distance travelled by a moving particle between successive collisions, the _mean free_ _path_  $T B D$ .

From kinetic theory, the mass flow rate of particles passing through a surface of unit area $T B D$ . The mean $T B D$ -velocity of particles crossing the plane from the $T B D$ -direction is $T B D$ ; similarly from the $T B D$ -direction, it is $T B D$ .

The net $T B D$ -momentum of the particles, positive on the $T B D$ side of the plane is then

$ T B D $

The net momentum is equivalent to the shear stress on the $T B D$ side of the plane, $T B D$ , as described in Sec. 2.6. By comparison with Eq. (6.16), the dynamic viscosity $T B D$ in terms of molecular properties. The kinematic viscosity is

$ T B D $
The original analysis of Maxwell#footnote[James Clerk Maxwell, _Illustrations of the dynamical theory of gases. Part I. On the motions and_ _collisions of perfectly elastic spheres_, 1860.] uses $T B D$ in Eq. (6.17). It was later recognised the average distance described by $T B D$ was larger due to _persistence of velocities_, _i.e._ a particle will sometimes maintain a path towards the plane after a collision.
A more thorough analysis#footnote[see Sydney Chapman and Thomas Cowling, _The mathematical theory of non-uniform gases (3rd ed.)_, 1970.] begins with the Boltzmann equation and applies the Chapman–Enskog expansion to first order in Knudsen number

$ T B D $
The analysis derives the conservation laws with Newtonian and Fourier constitutive models, _i.e._ Eq. (2.19), Eq. (2.51), Eq. (2.41) and Eq. (2.54). Treating particles as rigid spheres and collisions as elastic, yields a value $T B D$ , leaving a simple expression for viscosity which is
$ T B D $



== Turbulent mixing

Turbulent flow is characterised by significant mixing of fluid eddies as the Reynolds experiment in Sec. 6.1 shows. CFD simulations generally need to accommodate turbulent mixing since it influences the diffusion of mass, momentum and energy. While the fluid mixing by mass diffusion itself can be important, the effect on momentum diffusion is often critical because it impacts the calculation of viscous forces and, thus, the flow itself.

=== Eddy viscosity

Boussinesq was the first to devise a model for turbulence. He recognised the similarity between the random motion of both eddies in a turbulent fluid and particles at a molecular scale.

By analogy to kinetic theory in Sec. 6.10, shear stresses due to turbulence are caused by the net momentum, tangential to a plane, due to the motion of eddies. Boussinesq related this shear stress to the velocity gradient through an _eddy viscosity_ $T B D$ .#footnote[Joseph Boussinesq, _Essai sur la th__éorie des eaux courantes_, 1877.]

He presented the turbulent stress $T B D$ in tensor form, including a pressure. Kinetic theory relates pressure to fluctuations#footnote[Daniel Bernoulli, _Hydrodynamica, sive de viribus et motibus fluidorum commentarii_, 1738.] in particle velocity $T B D$ by $T B D$ ; the kinetic energy associated with the fluctuations is $T B D$ . Applying the same argument to velocity fluctuations $T B D$ due to turbulence, leads to a turbulent “pressure” $T B D$ , where $T B D$ is the _turbulent kinetic energy_ per unit mass.

By analogy with the Newtonian fluid model Eq. (2.41), the eddy viscosity model of Boussinesq, incorporating $T B D$ , is

$ T B D $
where $T B D$ is the viscous component of Reynolds stress and $T B D$ is the deformation rate tensor defined in Eq. (2.33). Inevitably Eq. (6.20) and Eq. (2.41) closely resemble one another.
The model of Eq. (6.20) requires some means of calculating $T B D$ . Kinetic theory gives a quantitative prediction of $T B D$ in Eq. (6.19) which led Boussinesq to hypothesise that $T B D$ , where $T B D$ and $T B D$ are a representative speed and length, respectively, with the speed $T B D$ relating to $T B D$ due to turbulence.

The turbulent viscosity can also be expressed as $T B D$ by absorbing the constant of proportionality within a characteristic speed $T B D$ and a _mixing length_ $T B D$ , discussed in Sec. 6.12.

Eddy viscosity and mixing length are useful concepts in turbulence modelling. However, it should be recognised that there are limitations in the analogy with kinetic theory, _e.g._:

- momentum is exchanged between submicroscopic particles through intermittent, discrete collisions, compared to the continuous interaction between eddies;
- the magnitude of random particle motions is generally equal in all directions, whereas the level of turbulent fluctuations can sometimes vary significantly with direction.




== Mixing length

By analogy with kinetic theory, the eddy viscosity $T B D$ can be expressed in terms of a characteristic speed $T B D$ and length $T B D$ . Prandtl produced the following model for the mixing length $T B D$ :#footnote[Ludwig Prandtl, _Bericht uber Untersuchungen zur ausgebildeten Turbulenz_, 1925.]

$ T B D $
The mixing length Eq. (6.21) is derived from the kinetic view of viscosity in Sec. 6.10. Turbulent fluctuations in the $T B D$ -direction, $T B D$ , carry fluid to the $T B D$ plane from a distance $T B D$ .
At $T B D$ , the turbulent fluctuations in the $T B D$ -direction, $T B D$ , correspond to the range of velocities at $T B D$ such that

$ T B D $
For a positive $T B D$ , a positive $T B D$ at $T B D$ corresponds to a negative $T B D$ , and vice versa. The Reynolds stress $T B D$ can then be constructed by defining a mixing length $T B D$ to replace $T B D$ , absorbing the constant of proportionality, to give
$ T B D $
From $T B D$ , we arrive at $T B D$ in Eq. (6.21).
The mixing length Eq. (6.21) is only effective as a turbulence model to calculate $T B D$ when the flow is simple enough that $T B D$ can be chosen appropriately.

Such an example is high $T B D$ , fully-developed flow through a pipe of radius $T B D$ . The mixing length $T B D$ , calculated from measured velocity profiles, follows a polynomial function of distance $T B D$ from the wall, given by#footnote[Johann Nikuradse, _Str__ömungsgesetze in rauhen Rohren_, 1933.]

$ T B D $

Notably, close to the wall, _e.g._ $T B D$ , the mixing length increases linearly according to

$ T B D $
where $T B D$ is Kármán’s constant.
At the centre of the pipe, $T B D$ where $T B D$ is pipe diameter. Estimates of $T B D$ are commonly cited for other simple examples, _e.g._ mixing layer, jet, flat plate boundary layer, _etc._, where $T B D$ is the characteristic length of the problem (radius in the case of a jet).



== Turbulent kinetic energy

The Reynolds stress $T B D$ in the ensemble-averaged momentum conservation Eq. (6.15) can be decomposed into viscous and pressure components. The turbulent pressure term, $T B D$ , is commonly subsumed within pressure $T B D$ , to give

$ T B D $
The viscous stresses can be combined assuming $T B D$ is modelled as a Newtonian fluid Eq. (2.41) and $T B D$ by the eddy viscosity model Eq. (6.20) to give
$ T B D $
where the _effective_ kinematic viscosity
$ T B D $
The effective viscosity represents momentum diffusion from the combined molecular and turbulent motions. The properties due to molecular motion are often described as _laminar_, _e.g._ laminar viscosity $T B D$ .
Averaging the momentum equation and introducing the eddy viscosity model creates one additional unknown $T B D$ . Additional models are required for $T B D$ to close the system of equations.

By considering $T B D$ from Sec. 6.11, the model for $T B D$ is typically decomposed into components representing velocity and length scales, $T B D$ and $T B D$ respectively.

The scale of $T B D$ corresponds to turbulent fluctuations $T B D$ , so it is reasonable to assume that $T B D$ . Since the field $T B D$ is representative of the $T B D$ component of $T B D$ , it is commonly adopted within turbulence models based on $T B D$ . Also, it is well captured by a suitable conservation equation.

=== Transport of turbulent kinetic energy

Conservation of turbulent kinetic energy $T B D$ can be written as:

$ T B D $
where the _turbulence generation_ $T B D$ is
$ T B D $
and $T B D$ is the effective diffusivity for $T B D$ . The equation is derived in a manner similar to Eq. (2.56) for conservation of specific internal energy $T B D$ , by (ensemble) averaging the separate energy contributions from $T B D$ and $T B D$ . While $T B D$ includes the kinetic energy of molecular motion, $T B D$ is the equivalent for eddy motion.
In Eq. (2.56), energy from bulk motion is passed to the submicroscopic scale as _heat_ by $T B D$ . In Eq. (6.28), it is converted to _turbulent energy_ by $T B D$ using Boussinesq’s $T B D$ from Eq. (6.20). The shear component $T B D$ provides the non-recoverable $T B D$ in Eq. (6.28) and the second ( $T B D$ ) term yields $T B D$ .

While $T B D$ transfers kinetic energy from the bulk flow to $T B D$ , the final term $T B D$ transfers $T B D$ on to $T B D$ as dissipated heat. Here, $T B D$ is the _turbulent dissipation rate_ per unit mass from in Sec. 6.6 which, from an ensemble-averaged derivation of Eq. (6.28), is

$ T B D $
Finally, diffusion of $T B D$ is represented by $T B D$ , where $T B D$ . This represents diffusion by both molecular and turbulent motions and interactions, including an adjustable coefficient $T B D$ which is usually set to 1.
The $T B D$ equation goes part of the way to closing our system of equations. However, it introduces an _additional unknown_, $T B D$ , and the model for $T B D$ still _requires length scale_ $T B D$ .



== Turbulent dissipation rate

A complete model for $T B D$ is still required to solve the ensemble-averaged momentum equation, _e.g._ Eq. (6.26). The discussion in Sec. 6.11 indicates $T B D$ is the product of $T B D$ and $T B D$ , requiring two models to represent each scale. Since $T B D$ , $T B D$ can represent the speed scale, modelled by Eq. (6.28).

The turbulent dissipation rate $T B D$ matches the rate of transfer of kinetic energy down the energy cascade $T B D$ , as discussed in Sec. 6.6. This applies to all turbulent scales including the larger mixing length scales, so $T B D$ . Substituting $T B D$ and $T B D$ expressions into $T B D$ yields turbulent viscosity

$ T B D $
where $T B D$ is a constant. From empirical data, $T B D$ , except within the viscous and buffer layers close to a wall, see Sec. 7.4.
=== Transport of turbulent dissipation rate

A model is now required for $T B D$ , both to calculate $T B D$ by Eq. (6.31) and to provide the remaining unknown in Eq. (6.28). The model can be provided by a transport equation for $T B D$

$ T B D $
where $T B D$ is the combined molecular and turbulent diffusion with an adjustable coefficient $T B D$ , usually set to 1.3. The remaining coefficients $T B D$ and $T B D$ are tuned to capture the behaviour of a range of flows.
The $T B D$ _-equation_, Eq. (6.32), can be derived in terms of statistical properties, replacing high-order terms in $T B D$ by models with coefficients $T B D$ , $T B D$ .#footnote[B. I. Davydov, _Statistical dynamics of an incompressible turbulent fluid_, _Dokl. Akad. Nauk SSSR_ *136*, 1961.]

Alternatively, it can be obtained by multiplying the principal variable, $T B D$ or $T B D$ , in each term of Eq. (6.28) by $T B D$ and introducing coefficients $T B D$ , $T B D$ and $T B D$ .

The $T B D$ term in Eq. (6.32) causes $T B D$ to increase with $T B D$ . This is logical since the generated turbulence moves down the energy cascade, so ultimately affects the rate of dissipation.

The $T B D$ term is justified by considering the _free decay of turbulence_. If the fluid stops moving ( $T B D$ ) and turbulence is no longer generated ( $T B D$ ), then (assuming constant $T B D$ and ignoring diffusion) Eq. (6.28) and Eq. (6.32) reduce to

$ T B D $
respectively. Integrating the combined equations yields $T B D$ where the “0” subscript indicates initial values.
Further integration yields a decay in $T B D$ over time to the power $T B D$ and a timescale $T B D$ , which is a reasonable approximation to real behaviour.



== Summary of turbulence

- Turbulence is a spatially complex distribution of vorticity which advects itself in a chaotic manner, Sec. 6.2.
- Vorticity is mostly generated at solid boundaries and evolves by advection, diffusion and stretching/compression, Sec. 6.3.
- Boundary layers are thus the main source of vorticity for turbulence, Sec. 6.4.
- Turbulence occurs when instabilities cause the vorticity to become chaotic, characterised by Reynolds Number, Sec. 6.1.


=== Scales of turbulence

- Length and time scales of turbulent eddies initially relate to the characteristic scales of the generated vortices, _e.g._ by shedding, Sec. 6.5.
- The eddies become progressively smaller until they reach a size where the dissipation of their kinetic energy is significant, Sec. 6.6.
- The process of eddy breakup involves a cascade of kinetic energy from larger to smaller eddies in which eddy structures become increasingly isotropic, Sec. 6.7.


=== Simulating turbulence

- The Kolmogorov microscales represent the smallest scales of turbulence at which eddies die out by viscous dissipation.
- They are surprisingly small, _e.g._ a factor of $T B D$ separates the length scales of the largest and smallest eddies.
- Usually, it is prohibitively costly to capture the smallest scales in a CFD simulation, Sec. 6.8.
- Large-eddy simulation captures the large eddies, while using modelling for the unresolved scales.
- Reynolds-averaged simulations provide faster solutions by solving equations for averaged field variables, Sec. 6.9.
- The averaged momentum equation includes the Reynolds stress which requires modelling.


=== Turbulence modelling concepts

- Momentum exchange at the molecular scale is characterised by viscosity and kinetic theory provides good quantitative predictions of viscosity for gases, Sec. 6.10.
- Comparing the random motion of eddies in a turbulent fluid with particles at a molecular scale, Boussinesq introduced the concept of eddy viscosity, Sec. 6.11.
- Eddy viscosity is not a fluid property but is proportional to a characteristic speed and length scale of the turbulence; models are required which represent each of these scales.
- The speed scale is represented by turbulent kinetic energy $T B D$ , described by a transport equation Eq. (6.28), Sec. 6.13.
- The $T B D$ -equation includes a term for its rate of dissipation $T B D$ ; a transport equation for $T B D$ , _e.g._ Eq. (6.32), provides a model for that term — which also represents the length scale of turbulence, Sec. 6.14.
- The concept of mixing length to characterise the length scale of turbulence provides additional modelling, particularly for the near-wall region, Sec. 6.12.

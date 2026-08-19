#import "../book.typ": *

#show: bookstyle.with(
  book-title: "Meshing for OpenFOAM",
  book-subtitle: "A guide for snappyHexMesh",
  authors: "Walter Dal'Maz Silva",
  references: "book.bib",
  margins: 1.75cm,
)

// ----------------------------------------------------------------------------
// PREAMBLE
// ----------------------------------------------------------------------------

#preamble[Preamble]

Geometry conception, cleaning, and meshing are crucial starting points for every successful simulation. One must put their best effort in this initial phase to ensure accurate downstream simulation results. It can add up to more than 80% of the total time spent in a new simulation case. This should be enough to convince anyone on spending more time improving these skills, and thus, read this guide.

Here we will focus on the use of _snappyHexMesh_, a component of the OpenFOAM framework for generating highly-complex meshes. This is a command line tool with a text-based interface, which can be intimidating to newcomers. We will show that although a large number of parameters are required for tuning the tool and conceiving a good quality mesh for simulation purposes, there is a special set of them that will make the major part of the job. Although _snappyHexMesh_ is a mesher of its own, it does not work alone, and its workflow generally starts from a background mesh generated with _blockMesh_, a tool that we will also need to learn the basics.

There are at least two major _flavors_ of OpenFOAM out there, OpenFOAM _Foundation_ and OpenFOAM _COM_ (from ESI Group). In many aspects they are the same, but in the last few years the gap has grown considerably, not only in the keywords, but also in the workflows themselves. This text is based on OpenFOAM Foundation version 13 (and we expect to keep it up-to-date with later versions). Be aware of this potential issue if you are constrained to use other flavors.

The methodology is progressive: we build upon simple cases, illustrating progressive improvements of the mesh. These may be quality improvements or general resolution improvements. You should keep in mind that _not all improvements_ must be used at all cases, at the risk of generating overly expensive simulations. One needs to understand what is the _Physics_ they are trying to capture in order to know when to stop _improving_ a mesh. When iterating over a case, we will assume you have read the preliminary steps for that case, so many aspects may be ommited to focus only in the new features.

We will discuss refinement controls, including the placement of _refinement regions_, both internal and external flow extraction, complex and moving geometries, as well as boundary layer resolution. @chapter-1 to @chapter-5 will present the fundamentals of each of the stages of meshing, while later chapters are exclusively guided examples. If this is not your first contact with the tool, you may skim the initial chapters and go straight into the examples. All examples are provided under _run/snappyHexMeshTutorials_.

// ----------------------------------------------------------------------------
// CONTENTS
// ----------------------------------------------------------------------------

#show: chapters
#include "_chapter1.typ"
#include "_chapter2.typ"
#include "_chapter3.typ"
#include "_chapter4.typ"
#include "_chapter5.typ"
#include "_chapter6.typ"

#show: appendix
#include "_appendix1.typ"

// ----------------------------------------------------------------------------
// EOF
// ----------------------------------------------------------------------------
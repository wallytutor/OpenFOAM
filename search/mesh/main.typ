// ----------------------------------------------------------------------------

#let book-title = "Meshing for OpenFOAM"
#let book-subtitle = "A guide for snappyHexMesh"
#let authors = "Walter Dal'Maz Silva"

// ----------------------------------------------------------------------------

#set document(title: book-title, author: authors)
#set text(lang: "en")

#set page(
  paper: "a5",
  margin: (
    inside: 2cm,
    outside: 2cm,
    top: 2cm,
    bottom: 2cm
  ),
  numbering: none,
  number-align: center,
)

#set par(justify: true, leading: 0.7em)

// Override the default "Section" used at all levels:
#show heading.where(level: 1): set heading(supplement: [Chapter])

// Ensure page break and vertical spacing before each chapter:
#show heading.where(level: 1): it => [
  #pagebreak(weak: true)
  #v(3em)
  #it
  #v(1.5em)
]

// Define preamble function to create unnumbered and unlisted heading:
#let preamble(title) = {
  heading(title, level: 1, numbering: none, outlined: false)
}

// ----------------------------------------------------------------------------
// COVER
// ----------------------------------------------------------------------------

#align(center + horizon)[
  #text(16pt, weight: "bold")[#book-title] \
  #v(1em)
  #text(16pt)[#book-subtitle] \
  #v(4em)
  #text(14pt)[#authors]
]

// ----------------------------------------------------------------------------
// TOC
// ----------------------------------------------------------------------------

#pagebreak()

#set page(numbering: "1")
#counter(page).update(1)

#outline(indent: 1.5em)

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

#set heading(numbering: "1.1")

#include "_chapter1.typ"
#include "_chapter2.typ"
#include "_chapter3.typ"
#include "_chapter4.typ"
#include "_chapter5.typ"

// ----------------------------------------------------------------------------
// LISTS
// ----------------------------------------------------------------------------

#bibliography(
  "references.bib",
  style: "ieee",
  title: "References"
)

#outline(
  title: [List of Figures],
  target: figure.where(kind: image),
)

#outline(
  title: [List of Tables],
  target: figure.where(kind: table),
)

// ----------------------------------------------------------------------------
// EOF
// ----------------------------------------------------------------------------
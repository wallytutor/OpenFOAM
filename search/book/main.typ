// ----------------------------------------------------------------------------

#let book-title = "Notes on Computational Fluid Dynamics"
#let book-subtitle = "General Principles Refactored"
#let authors = "Christopher Greenshields and Henry Weller"

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
  #v(1em)
  #text(14pt)[With notes by Walter Dal'Maz Silva]
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

Notes on Computational Fluid Dynamics (CFD) was written for people who use CFD in their work, research or study, providing essential knowledge to perform CFD analysis with confidence. It offers a modern perspective on CFD with the finite volume method, as implemented in OpenFOAM and other popular general-purpose CFD software. Fluid dynamics, turbulence modelling and boundary conditions are presented alongside the numerical methods and algorithms in a series of short, digestible topics, or notes, that contain complete, concise and relevant information. The book benefits from the experience of the authors: Henry Weller, core developer of OpenFOAM since writing its first lines in 1989; and, Chris Greenshields, who has delivered over 650 days of CFD training with OpenFOAM.

// ----------------------------------------------------------------------------
// CONTENTS
// ----------------------------------------------------------------------------

#set heading(numbering: none)

#include "_preface.typ"
#include "_symbols.typ"

#set heading(numbering: "1.1")

#include "_chapter1.typ"
#include "_chapter2.typ"
#include "_chapter3.typ"
#include "_chapter4.typ"
#include "_chapter5.typ"
#include "_chapter6.typ"
#include "_chapter7.typ"
#include "_chapter8.typ"

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
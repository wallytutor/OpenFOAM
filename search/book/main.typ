#import "../book.typ": bookstyle, preamble

#show: bookstyle.with(
  book-title: "Notes on Computational Fluid Dynamics",
  book-subtitle: "General Principles Refactored",
  authors: "Christopher Greenshields and Henry Weller",
  book-comment: "With notes by Walter Dal'Maz Silva",
  references: "book.bib",
  margins: 2.5cm,
)

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
// EOF
// ----------------------------------------------------------------------------
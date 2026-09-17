#import "@preview/glossarium:0.5.10": (
  make-glossary, register-glossary, print-glossary, gls
)

= Symbols

#show: make-glossary

#let user-print-glossary(entries, show-all: false) = {
  let user-print-title(entry) = text(weight: "regular", entry.short)
  let user-print-back-references(entry, deduplicate: false) = ""

  print-glossary(
    entries,
    user-print-title: user-print-title ,
    user-print-back-references: user-print-back-references,
    show-all: show-all,
  )
}

== General variables

#let general-variables = (
  (
    key: "x",
    short: $x$,
    description: [Cartesian axis direction $x$]
  ),
  (
    key: "y",
    short: $y$,
    description: [Cartesian axis direction $y$]
  ),
  (
    key: "z",
    short: $z$,
    description: [Cartesian axis direction $z$]
  ),
  (
    key: "r",
    short: $r$,
    description: [Spherical radial coordinate $r$]
  ),
  (
    key: "theta",
    short: $theta$,
    description: [Spherical polar coordinate $theta$]
  ),
  (
    key: "phi",
    short: $phi$,
    description: [Spherical azimual coordinate $phi$]
  ),
  (
    key: "V",
    short: $V$,
    description: [Volume or control volume $V$]
  ),
  (
    key: "S",
    short: $S$,
    description: [Surface or control surface $S$]
  ),
  // (
  //   key: "",
  //   short: $$,
  //   description: []
  // ),
)

#register-glossary(general-variables)
#user-print-glossary(general-variables)

== Geometry and mesh

== Numerical methods

== Physical quantites (with SI units)

#let symbol-entries = (
  (
    key: "pressure",
    short: $p$,
    description: [Pressure ($"Pa"$) or kinematic pressure ($"m"^2 "s"^(-2)$)],
  ),
  (
    key: "velocity",
    short: $bold(u)$,
    description: [Velocity vector ($"m" "s"^(-1)$)]

  )
)

#register-glossary(symbol-entries)
#user-print-glossary(symbol-entries)

== Relational symbols

// XXX using rel-\d was my strategy to enforce ordering, as sorting is done
// by key, not by short or description. There is no option to control this.
#let relational-symbols = (
  (
    key: "rel-1-equal",
    short: $a = b$,
    description: [$a$ and $b$ are equal]
  ),
  (
    key: "rel-2-approx",
    short: $a approx b$,
    description: [$a$ and $b$ are approximately equal]
  ),
  (
    key: "rel-3-almost",
    short: $a tilde.eq b$,
    description: [$a$ and $b$ are almost exactly equal]
  ),
  (
    key: "rel-4-equiv",
    short: $a equiv b$,
    description: [$a$ and $b$ are equivalent]
  ),
  (
    key: "rel-5-assigned",
    short: $a colon.eq b$,
    description: [$a$ is assigned the value of $b$]
  ),
  (
    key: "rel-6-propto",
    short: $a prop b$,
    description: [$a$ is proportional to $b$]
  ),
  (
    key: "rel-7-order",
    short: $a tilde b$,
    description: [$a$ is of the order of magnitude of $b$]
  )
)

#register-glossary(relational-symbols)
#user-print-glossary(relational-symbols, show-all: true)
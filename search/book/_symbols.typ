#import "@preview/glossarium:0.5.10": (
  make-glossary, register-glossary, print-glossary, gls
)

= Symbols

#show: make-glossary

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

#print-glossary(
  symbol-entries,
  user-print-title: entry => text(weight: "regular", entry.short),
  user-print-back-references: (entry, deduplicate: false) => "",
)

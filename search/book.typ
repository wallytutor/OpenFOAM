#let bookstyle(
  book-title: "",
  book-subtitle: "",
  book-comment: "",
  authors: "",
  margins: 1.75cm,
  paper: "a5",
  references: "",
  contents,
) = {
  set document(title: book-title, author: authors)
  set text(lang: "en")

  set page(
    paper: paper,
    margin: (
      inside:  margins,
      outside: margins,
      top:     margins,
      bottom:  margins
    ),
    numbering: none,
    number-align: center,
  )

  // Set default leading and justify text:
  set par(justify: true, leading: 0.7em)

  // Override the default "Section" used at all levels:
  show heading.where(level: 1): set heading(supplement: [Chapter])

  // Ensure page break and vertical spacing before each chapter:
  show heading.where(level: 1): it => [
    #pagebreak(weak: true)
    #v(3em)
    #it
    #v(1.5em)
  ]

  // Code styling
  show raw.where(block: true): it => {
    // Smaller font size:
    set text(size: 8pt)

    // Add line numbrers:
    show raw.line: line => {
      box(
        width: 2em,
        align(
          right,
          text(fill: luma(130), str(line.number))
        )
      )
      h(0.8em)
      line.body
    }

    // Add background and rounded corners:
    block(
      fill: rgb("#f6f8fa"),
      inset: 10pt,
      radius: 4pt,
      width: 100%,
      it,
    )
  }

  show figure.where(kind: raw): set align(left)
  show figure.where(kind: raw): set figure(supplement: [Listing])

  // Cover page
  align(center + horizon)[
    #text(16pt, weight: "bold")[#book-title] \
    #v(1em)
    #text(16pt)[#book-subtitle] \
    #v(4em)
    #text(14pt)[#authors]
    #v(2em)
    #text(12pt)[#book-comment]
  ]

  // Table of contents
  pagebreak()

  set page(numbering: "1")
  counter(page).update(1)

  // Contents
  outline(indent: 1.5em)
  contents

  // References
  if references != "" {
    bibliography(
      references,
      style: "ieee",
      title: "References"
    )
  }

  outline(
    title: [List of Figures],
    target: figure.where(kind: image),
  )

  outline(
    title: [List of Tables],
    target: figure.where(kind: table),
  )

  outline(
    title: [List of Listings],
    target: figure.where(kind: raw),
  )
}

#let preamble(title) = {
  heading(title, level: 1, numbering: none, outlined: false)
}
#let stylemain(margins: 1.75cm, paper: "a5", body) = {
  set page(
    paper: paper,
    margin: (
      inside: margins,
      outside: margins,
      top: margins,
      bottom: margins,
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
  body
}

#let stylecode(body) = {
  show raw.where(block: true): it => {
    // Smaller font size:
    set text(size: 8pt)

    // Add line numbrers:
    show raw.line: line => {
      box(
        width: 2em,
        align(
          right,
          text(fill: luma(130), str(line.number)),
        ),
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
  body
}

#let coverpage(
    book-title: "",
    book-subtitle: "",
    book-comment: "",
    authors: "",
  ) = {
  align(center + horizon)[
    #text(16pt, weight: "bold")[#book-title] \
    #v(1em)

    #if book-subtitle != "" {
      text(16pt)[#book-subtitle]
      v(4em)
    }

    #if authors != "" {
      text(14pt)[#authors]
      v(2em)
    }

    #if book-comment != "" {
      text(12pt)[#book-comment]
    }
  ]
}

#let lof() = {
  outline(
    title: [List of Figures],
    target: figure.where(kind: image),
  )
}

#let lot() = {
  outline(
    title: [List of Tables],
    target: figure.where(kind: table),
  )
}

#let lol() = {
  outline(
    title: [List of Listings],
    target: figure.where(kind: raw),
  )
}

#let enableoutlines() = {
  lof()
  lot()
  lol()
}

#let bookstyle(
    book-title: "",
    book-subtitle: "",
    book-comment: "",
    authors: "",
    margins: 1.75cm,
    paper: "a5",
    lang: "en",
    references: "",
    finallists: true,
    contents,
    ) = {
  set text(lang: lang)
  set document(title: book-title, author: authors)

  show: stylemain.with(margins: margins, paper: paper)
  show: stylecode

  coverpage(
    book-title: book-title,
    book-subtitle: book-subtitle,
    book-comment: book-comment,
    authors: authors,
  )

  pagebreak()
  counter(page).update(1)
  outline(indent: 1.5em)

  set page(numbering: "1")
  contents

  if references != "" {
    bibliography(
      references,
      style: "ieee",
      title: "References",
    )
  }
  if finallists {
    enableoutlines()
  }
}

#let preamble(title) = {
  heading(title, level: 1, numbering: none, outlined: false)
}

#let chapters(body) = {
  counter(heading).update(0)
  set heading(
    numbering: "1.1",
    supplement: [Chapter],
  )
  body
}

#let appendix(body) = {
  counter(heading).update(0)
  set heading(
    numbering: "A.1",
    supplement: [Annex],
  )
  body
}

// --------------------------------------------------------------------------
// Styles
// --------------------------------------------------------------------------

#let styleheading(body) = {
  // Override the default "Section" used at all levels:
  show heading.where(level: 1): set heading(supplement: [Chapter])

  // Part headings have no numbering (so they don't increment chapter counter)
  show heading.where(label: <part-heading>): set heading(numbering: none)

  // Handle headings with support to Part pages
  show heading: it => {
    if it.at("label", default: none) == <part-heading> {
      // Hide heading markup on the page so only the centered layout renders
      none
    } else {
      // Keep unmodified for other heading levels
      it
    }
  }

  // Ensure page break and vertical spacing before each chapter:
  show heading.where(level: 1): it => {
    if it.at("label", default: none) == <part-heading> {
      none
    } else [
      #pagebreak(weak: true)
      #v(3em)
      #it
      #v(1.5em)
    ]
  }

  body
}

#let stylefigure(body) = {
  // Apply a smaller font size for all captions; notice that we enforce
  // size unconditionally using an explicit content transformer:
  // show figure.caption: set text(size: 8.5pt) --> weak enforcement
  show figure.caption: it => text(size: 8.5pt)[
    #text(weight: "semibold", style: "normal")[
      #it.supplement
      #context it.counter.display(it.numbering):
    ]
    #it.body
  ]

  // Code listings specific adjustments:
  show figure.where(kind: raw): set align(left)
  show figure.where(kind: raw): set figure(supplement: [Listing])

  body
}

#let styleraw(body) = {
  show raw.where(block: true): it => {
    // Add line numbers:
    show raw.line: line => {
      box(
        width: 1.2em,
        align(
          right,
          text(fill: rgb("#7a9bb8"), str(line.number)),
        ),
      )
      h(0.8em)
      line.body
    }

    // Add background and rounded corners:
    block(
      fill: rgb("#f6f8fa"),
      inset: 5pt,
      radius: 0pt,
      width: 100%,

      // Size MUST be wrapped here
      text(size: 6pt, it),
    )
  }

  body
}

// --------------------------------------------------------------------------
// Cover page
// --------------------------------------------------------------------------

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

// --------------------------------------------------------------------------
// outline functions
// --------------------------------------------------------------------------

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

// --------------------------------------------------------------------------
// sectioning functions
// --------------------------------------------------------------------------

#let part-counter = counter("part")

#let preamble(title) = {
  heading(title, level: 1, numbering: none, outlined: false)
}

#let part(title) = {
  pagebreak(to: "odd")
  part-counter.step()

  // Create an outlined heading targeting parts
  [= #title <part-heading>]

  // Render the dedicated divider page
  context {
    let num = part-counter.display("I")
    align(center + horizon)[
      #text(size: 18pt, tracking: 2pt)[PART #num]
      #v(1.2em)
      #text(size: 26pt)[#title]
    ]
  }
  pagebreak()
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

// --------------------------------------------------------------------------
// bookstyle
// --------------------------------------------------------------------------

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

  show: styleheading
  show: stylefigure
  show: styleraw

  coverpage(
    book-title: book-title,
    book-subtitle: book-subtitle,
    book-comment: book-comment,
    authors: authors,
  )

  pagebreak()
  counter(page).update(1)

  show outline.entry: it => {
    if it.element.at("label", default: none) == <part-heading> {
      v(1.2em)
      let num = part-counter.at(it.element.location()).first()
      let roman = numbering("I", num)
      link(it.element.location())[
        #text(weight: "bold", size: 1.1em)[Part #roman: #it.element.body]
      ]
    } else {
      it
    }
  }

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

// --------------------------------------------------------------------------
// EOF
// --------------------------------------------------------------------------

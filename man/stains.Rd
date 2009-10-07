\name{stains}
\docType{data}
\alias{stains}
\title{Cytoband display information}
\description{
  This is a data set describing the display parameters used to plot
  cytoband data.
}
\usage{stains}
\format{A data frame with one row per cytoband type, and the following columns:
  \itemize{
    \item{type}{The cytoband type. This must correspond to the "stain"
      column in the cytoband data frame (see the \code{cytobands} documentation).}

    \item{bandcol}{The shade of gray used to colour the cytobands. A
      number between 0 (black) and 1 (white). Passed as the "col"
      argument to \code{rect}.}

    \item{textcol}{The shade of gray used for the cytoband text labels. A
      number between 0 (black) and 1 (white). Passed as the "col"
      argument to \code{text}.}

    \item{banddens}{The shading density to use for the band
    colour. Passed as the "density" argument to \code{rect}.}

    \item{bandbord}{The shade of gray used for the plotted cytoband
      borders. A number between 0 (black) and 1 (white). Passed as the
      "border" argument to \code{rect}.}
  }

}
\source{Developed based on the design of the idiogram Bioconductor package}
\keyword{datasets}

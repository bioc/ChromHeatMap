\name{cytobands}
\docType{data}
\alias{cytobands}
\title{Cytoband location information}
\description{
  This data set contains cytoband information for a range of species,
  taken directly from the UCSC genome annotation database. This data set
  is designed to be easily extendable to cover new species.
}
\usage{cytobands}
\format{A list of data frames, one per species, each with one row per
  cytoband and the following columns:

  \itemize{
    \item{chr}{The chromosome number for the cytoband, prefixed with 'chr'.}
    \item{start}{The start coordinate for the cytoband.}
    \item{end}{The end coordinate for the cytoband.}
    \item{band}{The cytoband number (i.e., the '23.3' in '1q23.3').}
    \item{stain}{The cytoband stain (see the \code{stains} data set).}
    \item{arm}{The chromosome arm for the cytoband (i.e., the 'q' in '1q23.3').}
  }

    The list names (i.e. \code{names(cytobands)}) should correspond to
    species names in the AnnotationDbi packages used.

}
\source{The UCSC genome annotation database: http://hgdownload.cse.ucsc.edu/downloads.html}
\keyword{datasets}

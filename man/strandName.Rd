\name{strandName}
\alias{strandName}

\title{Retrieve strand information from a ChrStrandMatrix object.}
\description{
This generic function simply returns the chromosome strand which name of all the chromosomes
represented by a given ChrStrandData object. Note that not every sample
associated with the object need have data from every chromosome.
}
\usage{
strandName(object)
}
\arguments{
  \item{object}{Object derived from class ChrStrandMatrix}
}
\value{
  \code{strandName(object)} returns the name of the strand from which
  the object data is taken.
}
\author{Tim F Rayner}

\seealso{\code{\link{ChrStrandMatrix-class}}}

\keyword{manip}

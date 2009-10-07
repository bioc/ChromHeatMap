\name{chrNames}
\alias{chrNames}

\title{Retrieve chromosome names from an object.}
\description{

  This generic function simply returns the names of all the chromosomes
  represented by a given ChrStrandData or ChrStrandMatrix object. Note
  that not every sample associated with a ChrStrandData object need have
  data from every chromosome.

}
\usage{
chrNames(object)
}
\arguments{
  \item{object}{Object derived from class ChrStrandData or ChrStrandMatrix}
}
\value{
  \code{chrNames(object)} returns a character vector listing the chromosomes.
}
\author{Tim F Rayner}

\seealso{\code{\link{ChrStrandData-class}}}

\keyword{manip}

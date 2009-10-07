\name{makeChrStrandData}
\alias{makeChrStrandData}

\title{ Map a data matrix onto chromosome coordinates }

\description{

  Given an ExpressionSet, or a data matrix with row names corresponding
  to the probe IDs in an accompanying annotation package, this function
  returns a data structure that can be used with the \code{plotChrMap}
  function. This code is based on the Makesense method from the
  geneplotter package, extended to use both the CHRLOC and CHRLOCEND
  annotation environments from recent AnnotationDbi packages.

}
\usage{
makeChrStrandData(expr, lib)
}
\arguments{
  \item{expr}{The ExpressionSet or data matrix to remap.}
  \item{lib}{The name of the annotation package to use.}
}
\value{

  A ChrStrandData object suitable for use with plotChrMap.

}
\references{ geneplotter, annotate and AnnotationDbi packages }
\author{ Tim F Rayner }
\seealso{
  \code{\link{plotChrMap}}, \code{\link{ChrStrandData-class}}
}
\examples{
data('demo')
chrdata <- makeChrStrandData(exprs(ALLs.chr22), lib = "hgu95av2.db")
}
\keyword{array}

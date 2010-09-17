\name{grabChrMapProbes}
\alias{grabChrMapProbes}

\title{ Identify the probes or genes plotted using plotChrMap }
\description{

  Allows the user to interactively select regions of the plotChrMap
  heatmap, identifying all the probes or genes plotted in those regions.

}
\usage{
grabChrMapProbes( plotmap )
}
\arguments{
  \item{plotmap}{The output of the \code{plotChrMap} function.}
}
\details{

  This function takes the output of the \code{plotChrMap} function and
  uses it to identify the probes or genes responsible for the signals plotted on
  the \code{plotChrMap} heatmap. It asks the user to select two points
  on either side of the heatmap bands of interest (specifically,
  boundary for inclusion of a given band is its left-hand edge), and
  returns a vector of probe/gene identifiers. This can be passed directly to
  AnnotationDbi::mget to yield gene symbols and other annotation.

  Note that the plotting area layout() and par() values are not reset on
  exit, so that this function can be reused as many times as is desired.

}
\value{

  A character vector of probe/gene identifiers. If multiple identifiers
  have been averaged into a single band these identifiers will be string
  concatenated, separated by semicolons. The start, end and interval
  arguments to \code{plotChrMap} can be used in such cases to plot the
  data at a higher resolution, splitting such loci into separate bands.

}
\author{ Tim F Rayner }
\seealso{
  \code{\link{plotChrMap}}
}
\examples{
data('demo')
plotmap <- plotChrMap(chrdata, '22', cytoband='q11.23')
probes <- grabChrMapProbes(plotmap)
library('hgu95av2.db')
genes <- mget(probes, hgu95av2SYMBOL, ifnotfound=NA)
}
\keyword{hplot}


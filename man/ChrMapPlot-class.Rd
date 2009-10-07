\name{ChrMapPlot}
\docType{class}
\alias{class:ChrMapPlot}
\alias{ChrMapPlot}
\alias{ChrMapPlot-class}

\alias{show,ChrMapPlot-method}
\alias{summary,ChrMapPlot-method}

\title{
  Class containing a mapping between plot location and probe identifier.
}

\description{
  ChrMapPlot objects are generated as an output from the main
  \code{plotChrMap} function, which users can then pass to the
  \code{grabChrMapProbes} function.
}
\section{Creating Objects}{

  Objects of this class are created using the \code{plotChrMap}
  function:

  \code{plotChrMap(chrdata, '22')}
}
\section{Slots}{

   \describe{

     \item{labels}{An array of probe identifiers, with names corresponding
       to chromosome coordinates}.

     \item{start}{The leftmost interval number (most usually 1).}

     \item{end}{The rightmost interval number.}
  }
}
\section{Methods}{

  Standard generic methods:
  \describe{
    \item{\code{show(ChrMapPlot)}}{Generates a short description of
      the ChrMapPlot object.}
  }
}

\author{Tim F Rayner}

\seealso{
  \code{\link{plotChrMap}}, \code{\link{grabChrMapProbes}}.
}

\examples{
data('demo')
plotmap <- plotChrMap(chrdata, '22', cytoband='q11.23')
probes <- grabChrMapProbes(plotmap)
library('hgu95av2.db')
genes <- mget(probes, hgu95av2SYMBOL, ifnotfound=NA)
}

\keyword{classes}

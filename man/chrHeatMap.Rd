\name{chrHeatMap}
\alias{chrHeatMap}

\title{ Plot ChrStrandMatrix objects as heat maps along a chromosome }

\description{

  Plots a either one or two ChrStrandMatrix objects (typically
  constructed using the \code{createChrMatrix} function) as heat maps
  along a specified chromosome, optionally clustering samples and
  including an idiogram.

}
\usage{
chrHeatMap (strand.data, cytopaint.func=NULL, col = "heat.colors",
            start, end, breaks, RowSideColors, title=TRUE,
            margins = c(6, 6), cexCyto = 0.8, srtCyto=90, lmat = NULL, lhei = NULL,
            lwid = NULL, ...)
}
\arguments{

  \item{strand.data}{A ChrStrandMatrix object, or a list of such
    objects, one per strand to be plotted (or a single matrix for `both'
    strands), created using the \code{createChrMatrix} function.}

  \item{cytopaint.func}{A function closure taking a single argument, `boxwidth',
    and plotting its enclosed idiogram data at that width. See
    \code{plotChrMap} for the code used to generate this closure.}

  \item{col}{A vector of colors to use for the heat map, or the name of
    a function generating such a vector.}

  \item{start}{The starting genome coordinate for the plot.}

  \item{end}{The ending genome coordinate for the plot.}

  \item{breaks}{A vector of numeric break points indicating the boundaries
    between the \code{col} colors.}

  \item{RowSideColors}{A vector of colors to use for a color band
    indicating e.g. sample categories.}

  \item{title}{If TRUE, this causes the function to include default heat
    map subtitles indicating which chromosome and strand has been
    plotted. If FALSE or NULL, subtitles will left blank. If this
    argument is set to a character vector of the same length and order
    as strand.data its contents will be used as heat map subtitles.}

  \item{margins}{A numeric vector indicating the c(bottom, left) margins
    of the plot containing X and Y axes labels.}

  \item{cexCyto}{A positive number used to control the font size for the
    idiogram plot. For plots spanning just a few cytobands it may be
    worth setting this to a larger number, and \code{srtCyto}, below, to
    zero.}

  \item{srtCyto}{A number indicating the degree to which the idiogram
    text labels should be rotated. This defaults to 90 degrees, but for
    more detailed plots a setting of zero here often looks better.}

  \item{lmat}{An optional matrix to be passed to \code{layout}.}

  \item{lhei}{An optional vector of \code{layout} row heights.}

  \item{lwid}{An optional vector of \code{layout} row widths.}

  \item{...}{Additional arguments are passed to the \code{drawMapDendro} function.}

}
\details{

  Typically this function should not be called directly, but rather via
  the wrapper \code{plotChrMap} function. This function uses cytoband
  data from the UCSC genome annotation database and code adapted from
  the quantsmooth package to draw an idiogram of the chromosome, or a 
  subset thereof.

}
\value{
  This function is executed for its side effects.
}
\references{ lodplot and quantsmooth packages }
\author{ Tim F Rayner }
\seealso{
  \code{\link{plotChrMap}}, \code{\link{createChrMatrix}}, \code{\link{drawMapDendro}}
}

\examples{
data('demo')
stranddata <- createChrMatrix( chrdata, chr=22, strand='forward', start=21925000, end=24300000 )
chrHeatMap(stranddata)
}
\keyword{hplot}

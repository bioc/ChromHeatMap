\name{plotChrMap}
\alias{plotChrMap}

\title{ Plot data as an annotated heat map along a chromosome }
\description{

  Given a ChrStrandData object (produced by the \code{makeChrStrandData}
  function), this function plots a heat map of its data values along a
  specified chromosome, optionally clustering samples and including an
  idiogram.

}
\usage{
plotChrMap( data, chr, start = 1, end, subset = NULL,
            cytoband, interval = ceiling((end-start)/500),
            strands = c('forward', 'reverse'), ... )
}
\arguments{

  \item{data}{A ChrStrandData object, output from the \code{makeChrStrandData} function.}

  \item{chr}{Chromosomal id, chromosome to plot 1:22,X,Y.}

  \item{start}{Optional start chromosome position from which to commence plotting.}

  \item{end}{Optional end chromosome position.}

  \item{subset}{Optional numeric vector listing the samples from \code{data} to plot.}

  \item{cytoband}{Optional cytological band to plot (e.g. `q23').}

  \item{interval}{An optional interval size controlling the plot detail
    level.}

  \item{strands}{The chromosome strands to plot (a one- or two-element
    character vector, values `forward', `reverse', or `both').}

  \item{...}{Additional arguments are passed to the \code{chrHeatMap} function.}

}
\details{

  This function is used to plot ChrStrandData objects (the output of the
  \code{makeChrStrandData} function) as heatmaps arranged along genome
  coordinates. The default heat map will plot the entire forward strand
  for the chosen chromosome at the top of the figure, with an idiogram
  and the reverse strand below it. To plot both strands overlaid, use
  the \code{strands='both'} argument. Probe signals are averaged over a
  window size controlled by \code{interval}, such that the default
  length of each heat map segment is 1/500 the total heat map
  width. This can be varied as required to control the resolution of the
  plot. This function uses both the start and end chromosomal locations
  for each gene to plot heatmap positions, and as such will not work
  with older AnnotationDbi packages.

  See the related functions from this package for further plotting
  arguments which may be passed to this function. In particular, see the
  \code{drawMapDendro} documentation for arguments used to control
  sample clustering and plot axis font sizes, and \code{chrHeatMap} for
  arguments relating to the idiogram plot. Note that the plotting area
  layout() and par() values are not reset on exit, so that
  \code{grabChrMapProbes} can be subsequently used on the output.

  Idiogram plotting is currently only supported for data mapping to
  human, mouse and rat genomes. In principle this is extendable to any
  organism for which the UCSC genome browser includes cytoband
  information. Please contact the maintainer of this package for help in
  such cases.

}
\value{

  A ChrMapPlot object containing a list of probe identifiers mapped to
  their corresponding display locations, for use with
  \code{grabChrMapProbes}.
  
}
\references{ annotate package }
\author{ Tim F Rayner }
\seealso{
  \code{\link{drawMapDendro}}, \code{\link{chrHeatMap}}, \code{\link{makeChrStrandData}}, \code{\link{grabChrMapProbes}}
}
\examples{
data('demo')
plotChrMap(chrdata, '22', cytoband='q11', labRow=ALLs.chr22$mol.biol,
cexCol=0.8, cexCyto=1.2, srtCyto=0)
}
\keyword{hplot}


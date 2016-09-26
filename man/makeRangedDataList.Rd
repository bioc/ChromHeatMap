\name{makeRangedDataList}
\alias{makeRangedDataList}

\title{ Plot expression data as tracks in the UCSC genome browser }
\description{

  Creates a GRangesList object suitable for uploading to the UCSC
  genome browser using the rtracklayer package. 

}
\usage{
makeRangedDataList( data, chr, start = 1, end, genome, subset = NULL,
                    cytoband, plot=FALSE, session )
}
\arguments{

  \item{data}{A ChrStrandData object, output from the \code{makeChrStrandData} function.}

  \item{chr}{Chromosomal id, chromosome to plot 1:22,X,Y.}

  \item{start}{Optional start chromosome position from which to commence plotting.}

  \item{end}{Optional end chromosome position.}

  \item{genome}{The name of the genome from which the data coordinates
    are taken (e.g. "hg18"). Passed to \code{GenomicData} in the
    rtracklayer package.}

  \item{subset}{Optional numeric vector listing the samples from \code{data} to plot.}

  \item{cytoband}{Optional cytological band to plot (e.g. `q23').}

  \item{plot}{An optional flag indicating whether to automatically plot
    the resulting GRangesList on the UCSC browser or not.}

  \item{session}{An optional rtracklayer UCSCSession object. Ignored
    unless plot=TRUE.}
}
\details{

  This function is used to create GRangesList objects from ChrStrandData
  objects (see the \code{makeChrStrandData} function). If the \code{plot}
  argument is set to TRUE, the data is also uploaded to a UCSC browser
  session using default settings. See the rtracklayer package for more
  information on RangedData and UCSCSession objects.

}
\value{

  A GRangesList object containing the data for the specified genome
  region. See the rtracklayer package for more information on this
  object class.

}
\references{ rtracklayer package }
\author{ Tim F Rayner }
\seealso{
  \code{\link{makeChrStrandData}},
  \code{\link[GenomicRanges:GRangesList-class]{GRangesList}}
  \code{\link{plotChrMap}},
}
\examples{
data('demo')
r <- makeRangedDataList( data=chrdata, chr=22, cytoband='q11.23', genome='hg18' )
}

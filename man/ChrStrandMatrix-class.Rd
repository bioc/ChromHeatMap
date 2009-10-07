\name{ChrStrandMatrix}
\docType{class}
\alias{class:ChrStrandMatrix}
\alias{ChrStrandMatrix}
\alias{ChrStrandMatrix-class}

\alias{show,ChrStrandMatrix-method}
\alias{summary,ChrStrandMatrix-method}
\alias{chrNames,ChrStrandMatrix-method}
\alias{strandName,ChrStrandMatrix-method}
\alias{sampleNames,ChrStrandMatrix-method}
\alias{featureNames,ChrStrandMatrix-method}
\alias{exprs,ChrStrandMatrix-method}

\title{
  Class to contain data associated with genome locations for a specific chromosome.
}

\description{
  Container for chromosome-specific subsets of data selected from an
  genome-wide ChrStrandData object, suitable for use with
  \code{chrHeatMap}.
}

\section{Creating Objects}{

  Typically, objects of this class are created and used internally by
  the \code{createChrMatrix} and \code{chrHeatMap} functions. Objects
  can be created in a similar fashion by end-users:

  \code{createChrMatrix(chrdata, chr=22, strand='forward', start=21925000,
    end=24300000, interval=5000)}

  Note that this function may combine data from multiple probes (taking
  the mean) into a single chromosomal locus based on the size of the
  specified interval. If this happens the combined probe identifiers are
  concatenated in the output object, separated by a semicolon.
  
}

\section{Slots}{

   \describe{

     \item{data}{The data matrix, arranged with samples in columns and
       genomic locations in rows.}

     \item{probeID}{An array of probe identifiers associated with the
       data. The names attached to this array correspond with chromosome
       coordinate (specifically, the starting coordinates, i.e. the
       left-hand edges). These identifiers will ultimately be returned by e.g. the
       \code{grabChrMapProbes} function.}

     \item{chr}{The chromosome name or number.}

     \item{strand}{The chromosome strand ('forward', 'reverse' or 'both').}

     \item{start}{The starting chromosome coordinates for each genomic location.}

     \item{end}{The ending chromosome coordinates for each genomic
       location.}
  }
}
\section{Methods}{

  Class-specific methods.
  \describe{

    \item{\code{chrNames(ChrStrandMatrix)}}{Returns the name of the
      chromosome for the object.}

    \item{\code{strandName(ChrStrandMatrix)}}{Returns the chromosome
      strand for the object.}

    \item{\code{sampleNames(ChrStrandMatrix)}}{Returns the names of the
      samples associated with the object.}    

    \item{\code{featureNames(ChrStrandMatrix)}}{Returns the probe
      identifiers associated with the object.}    

    \item{\code{exprs(ChrStrandMatrix)}}{Returns the chromosome-specific
      data matrix for the object.}    
  }

  Standard generic methods:
  \describe{

    \item{\code{show(ChrStrandMatrix)}}{Generates a short description of
      the ChrStrandMatrix object.}

    \item{\code{summary(ChrStrandMatrix)}}{Generates a summary of the data
      available for each sample in the ChrStrandMatrix object.}
  }
}

\author{Tim F Rayner}

\seealso{
  \code{\link{createChrMatrix}}, \code{\link{ChrStrandData-class}}.
}

\examples{
data('demo')
stranddata <- createChrMatrix( chrdata, chr=22, strand='forward', start=21925000, end=24300000 )
}

\keyword{classes}

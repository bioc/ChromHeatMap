\name{ChrStrandData}
\docType{class}
\alias{class:ChrStrandData}
\alias{ChrStrandData}
\alias{ChrStrandData-class}

\alias{show,ChrStrandData-method}
\alias{summary,ChrStrandData-method}
\alias{annotation,ChrStrandData-method}
\alias{chrNames,ChrStrandData-method}
\alias{sampleNames,ChrStrandData-method}

\title{
  Class to contain data associated with chromosome coordinates across a
  whole genome.
}

\description{
  Container for data from high-throughput assays mapped to chromosome locations.
}
\section{Creating Objects}{

  The most convenient way to create a ChrStrandData object is to use the
  \code{makeChrStrandData} function, which can be used to convert data
  stored in either an \code{ExpressionSet} or data frame into a
  \code{ChrStrandData} object:

  \code{makeChrStrandData(ALL, lib = "hgu95av2.db")}
}
\section{Slots}{

   \describe{

     \item{data}{a 'list', whose components correspond to samples in the
       same order as appearing in the columns of 'expr'. Each component
       is also a 'list', named by chromosomes "1"-"22", "X" and
       "Y". Each named component is again a 'list' with two elements
       named "posS" and "negS", corresponding to the forward and
       reverse strands of a chromosome, each of which is a list
       containing start coordinates ("x"), end coordinates("xe") and the
       corresponding data values ("y").}

     \item{lib}{A string giving the name of the annotation data package
       to use.}

     \item{chrs}{The list of chromosomes represented in the object.}
  }
}
\section{Methods}{

  Class-specific methods.
  \describe{

    \item{\code{annotation(ChrStrandData)}}{Returns the name of the
      AnnotationDbi library used to annotate the object.}

    \item{\code{chrNames(ChrStrandData)}}{Returns a list of the 
      chromosomes represented in the object.}

    \item{\code{sampleNames(ChrStrandData)}}{Returns the names of the
      samples associated with the object.}    
  }
  
  Standard generic methods:
  \describe{

    \item{\code{show(ChrStrandData)}}{Generates a short description of
      the ChrStrandData object.}

    \item{\code{summary(ChrStrandData)}}{Generates a summary of the data
      available for each chromosome in the ChrStrandData object.}
  }
}

\author{Tim F Rayner}

\seealso{
  \code{\link{makeChrStrandData}}, \code{\link{ChrStrandMatrix-class}}.
}

\examples{
data('demo')
chrdata <- makeChrStrandData(exprs(ALLs.chr22), lib = "hgu95av2.db")
}

\keyword{classes}

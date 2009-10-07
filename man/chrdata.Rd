\name{chrdata}
\docType{data}
\alias{chrdata}
\title{The ALLs.chr22 ExpressionSet, reformatted as a ChrStrandData object}
\description{
  This is a greatly reduced subset of the Chiaretti et al. ALL data
  set (available in its entirety as the Bioconductor ALL package). The
  data in this subset consist of microarrays from 15 different
  individuals with acute lymphoblastic leukemia (ALL). The data are
  further restricted to chromosome 22 only. This data set is intended
  for demonstration purposes only. See the documentation for
  \code{makeChrStrandData} for a description of the ChrStrandData object
  format. This format directly associates the ExpressionSet data with
  chromosome location, speeding up retrieval of data during heat map plotting.
}
\usage{chrdata}
\format{A ChrStrandData object}
\source{The ALL Bioconductor data package}
\references{
  Sabina Chiaretti, Xiaochun Li, Robert Gentleman, Antonella Vitale,
  Marco Vignetti, Franco Mandelli, Jerome Ritz, and Robin Foa Gene
  expression profile of adult T-cell acute lymphocytic leukemia
  identifies distinct subsets of patients with different response to
  therapy and survival. Blood, 1 April 2004, Vol. 103, No. 7.
}
\keyword{datasets}

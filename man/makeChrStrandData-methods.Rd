\name{makeChrStrandData-methods}
\docType{methods}
\alias{makeChrStrandData-methods}
\alias{makeChrStrandData,ExpressionSet-method}
\alias{makeChrStrandData,matrix-method}
\title{Map a data matrix onto chromosome coordinates}
\description{

  Given a data matrix with row names corresponding to the probe or gene IDs in
  an accompanying annotation package, returns a data structure that can
  be used with the \code{plotChrMap} function. Based on the Makesense
  method from the geneplotter package.

}
\section{Methods}{
\describe{

\item{expr = "ExpressionSet"}{ Given an ExpressionSet object, returns a
  ChrStrandData object. }

\item{expr = "matrix"}{ Given a matrix object (where
  \code{rownames(expr)} yields the probe or gene identifiers used by the
  annotation package), returns a ChrStrandData object. }

}}
\keyword{methods}

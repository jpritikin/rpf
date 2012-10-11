\docType{package}
\name{An introduction}
\alias{An introduction}
\alias{An introduction-package}
\title{rpf - Response Probability Functions}
\description{
  The purpose of this package is to factor out logic and
  math common to Item Response Theory fitting, diagnostics,
  and analysis.  It is envisioned as core support code
  suitable for more specialized IRT packages to build upon.
}
\details{
  This package provides optimized, low-level functions to
  map parameters to response probabilities for dichotomous
  (1PL, 2PL and 3PL) \code{\link{rpf.drm}} and polytomous
  (graded response \code{\link{rpf.grm}}, partial
  credit/generalized partial credit \code{\link{rpf.gpcm}},
  nominal \code{\link{rpf.nrm}}, and multiple-choice model
  \code{\link{rpf.mcm}}) items. Both unidimensional and
  multidimensional versions of the models are available.

  Item model parameters are passed around as a numeric
  vector. A 1D matrix is also acceptable. Regardless of
  model, parameters are always ordered as follows:
  discrimination ("a"), difficulty ("b"), and guessing
  ("c").

  This package could also accrete functions to support
  plotting (but not the actual plot functions).
}
\section{Warning}{
  The API is not stable at this time. You have been warned.
}
\seealso{
  See \code{\link{rpf.rparam}} and
  \code{\link{rpf.startingParam}} to create item
  parameters.
}


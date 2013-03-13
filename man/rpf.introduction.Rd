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
  ("c"). If person ability ranges from low negative to high
  positive then probabilities are output from incorrect to
  correct. That is, a low ability person (e.g., ability =
  -2) will be more likely to get an item incorrect than
  correct. For example, a dichotomous model that returns
  [.25, .75] indicates a probability of .25 for incorrect
  and .75 for correct.  A polytomous model will have the
  most incorrect probability at index 1 and the most
  correct probability at the maximum index.

  All models are always in the logistic metric. To obtain
  normal ogive discrimination parameters, divide slope
  parameters by \code{\link{rpf.ogive}}. Item models are
  estimated in slope-intercept form unless the traditional
  parameterization is specifically requested.

  This package could also accrete functions to support
  plotting (but not the actual plot functions).
}
\section{Warning}{
  The API is not stable at this time. You have been warned.
  In particular, I anticipating transposing some of the
  input and output matrices.
}
\references{
  Thissen, D. and Steinberg, L. (1986). A taxonomy of item
  response models. \emph{Psychometrika 51}(4), 567-577.
}
\seealso{
  See \code{\link{rpf.rparam}} to create item parameters.
}


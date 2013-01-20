\name{rpf.gpcm}
\alias{rpf.gpcm}
\title{Create a generalized partial credit model and associated hyperparameters.}
\usage{
  rpf.gpcm(numOutcomes = 2, dimensions = 1,
    multidimensional)
}
\arguments{
  \item{numOutcomes}{The number of choices available}

  \item{dimensions}{the number of dimensions}

  \item{multidimensional}{whether to use a multidimensional
  model. Defaults to \code{TRUE} when \code{dimensions>1}
  and \code{FALSE} when \code{dimensions==1}.}
}
\value{
  an item model
}
\description{
  Create a generalized partial credit model and associated
  hyperparameters.
}
\references{
  Baker & Kim (2004). \emph{Item Response Theory: Parameter
  Estimation Techniques.} Marcel Dekker, Inc.

  Masters, G. N. (1982). A rasch model for partial credit
  scoring. \emph{Psychometrika, 47}(2), 149-174.

  Muraki, Eiji. (1992). A generalized partial credit model:
  Application of an EM algorithm. \emph{Applied
  Psychological Measurement, 16}, 159-176.
  doi:10.1177/014662169201600206
}


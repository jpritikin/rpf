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
  Bayesian priors are only used to generate plausible
  random parameters.
}
\references{
  Baker & Kim (2004). Item Response Theory: Parameter
  Estimation Techniques. Marcel Dekker, Inc.
}


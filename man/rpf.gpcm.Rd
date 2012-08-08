\name{rpf.gpcm}
\alias{rpf.gpcm}
\title{Create a generalized partial credit model and associated hyperparameters.}
\usage{
  rpf.gpcm(numOutcomes = 2, dimensions = 1, D = 1,
    multidimensional)
}
\arguments{
  \item{numOutcomes}{The number of choices available}

  \item{dimensions}{the number of dimensions}

  \item{D}{defaults to 1 or pass in
  \code{\link{rpf.ogive}}}

  \item{multidimensional}{whether to use a multidimensional
  model. Defaults to \code{TRUE} when \code{dimensions>1}
  and \code{FALSE} when \code{dimensions==1}.}
}
\value{
  an item model
}
\description{
  This function instantiates a generalized partial credit
  model. The discrimination prior defaults to the lognormal
  distribution with \code{meanlog=0} and \code{sdlog=.5}.
}
\details{
  It is not yet possible to further customize the Bayesian
  priors. The API will change before the 1.0 release.
}
\references{
  Baker & Kim (2004). Item Response Theory: Parameter
  Estimation Techniques. Marcel Dekker, Inc.
}


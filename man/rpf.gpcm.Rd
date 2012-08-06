\name{rpf.gpcm}
\alias{rpf.gpcm}
\title{Create a generalized partial credit model and associated hyperparameters.}
\usage{
  rpf.gpcm(numOutcomes = 2, D = 1)
}
\arguments{
  \item{numOutcomes}{The number of choices available}

  \item{D}{defaults to 1 or pass in the
  \code{\link{rpf.ogive}}}
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


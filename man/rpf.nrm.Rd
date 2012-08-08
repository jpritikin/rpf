\name{rpf.nrm}
\alias{rpf.nrm}
\title{Create a nominal response model and associated hyperparameters.}
\usage{
  rpf.nrm(numOutcomes = 2, dimensions = 1)
}
\arguments{
  \item{numOutcomes}{The number of choices available}

  \item{dimensions}{the number of dimensions}
}
\value{
  an item model
}
\description{
  This function instantiates a nominal response model. The
  discrimination prior defaults to the lognormal
  distribution with \code{meanlog=0} and \code{sdlog=.5}.
}
\details{
  It is not yet possible to further customize the Bayesian
  priors. The API will change before the 1.0 release.
}
\author{
  Jonathan Weeks <weeksjp@gmail.com>
}


\name{rpf.mcm}
\alias{rpf.mcm}
\title{Create a multiple-choice response model and associated hyperparameters.}
\usage{
  rpf.mcm(numOutcomes = 2, numChoices = 5, dimensions = 1)
}
\arguments{
  \item{numOutcomes}{the number of possible outcomes}

  \item{numChoices}{the number of choices available}

  \item{dimensions}{the number of dimensions}
}
\value{
  an item model
}
\description{
  This function instantiates a multiple-choice response
  model. The discrimination prior defaults to the lognormal
  distribution with \code{meanlog=0} and \code{sdlog=.5}.
  The guessing prior is the beta distribution. See the
  source code for details.
}
\details{
  It is not yet possible to further customize the Bayesian
  priors. The API will change before the 1.0 release.
}
\author{
  Jonathan Weeks <weeksjp@gmail.com>
}


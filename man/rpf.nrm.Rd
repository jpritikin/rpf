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
  This function instantiates a nominal response model.
  Bayesian priors are only used to generate plausible
  random parameters.
}
\author{
  Jonathan Weeks <weeksjp@gmail.com>
}


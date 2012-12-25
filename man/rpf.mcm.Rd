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
  model. Bayesian priors are only used to generate
  plausible random parameters.
}
\author{
  Jonathan Weeks <weeksjp@gmail.com>
}


\name{rpf.grm}
\alias{rpf.grm}
\title{Create a graded response model and associated hyperparameters.}
\usage{
  rpf.grm(numOutcomes = 2, dimensions = 1,
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
  This function instantiates a graded response model.
  Bayesian priors are only used to generate plausible
  random parameters.
}
\details{
  The graded response model was designed for a item with a
  series of dependent parts where a higher score implies
  that easier parts of the item were surmounted. If there
  is any chance your polytomous item has independent parts
  then consider \code{\link{rpf.gpcm}}.
}
\author{
  Jonathan Weeks <weeksjp@gmail.com>
}


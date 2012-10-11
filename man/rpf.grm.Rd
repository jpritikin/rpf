\name{rpf.grm}
\alias{rpf.grm}
\title{Create a graded response model and associated hyperparameters.}
\usage{
  rpf.grm(numOutcomes = 2, dimensions = 1, D = 1,
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
  This function instantiates a graded response model. The
  discrimination prior defaults to the lognormal
  distribution with \code{meanlog=0} and \code{sdlog=.5}.
}
\details{
  It is not yet possible to further customize the Bayesian
  priors. The API will change before the 1.0 release.

  The graded response model was designed for a item with a
  series of dependent parts where a higher score implies
  that easier parts of the item were surmounted. If your
  polytomous item has a number of independent parts then
  consider \code{\link{rpf.gpcm}}.
}
\author{
  Jonathan Weeks <weeksjp@gmail.com>
}


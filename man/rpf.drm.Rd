\name{rpf.drm}
\alias{rpf.drm}
\title{Create a dichotomous response model and associated hyperparameters.}
\usage{
  rpf.drm(D = 1, numChoices = 5)
}
\arguments{
  \item{D}{defaults to 1 or pass in the
  \code{\link{rpf.ogive}}}

  \item{numChoices}{the number of alternatives in the
  question}
}
\value{
  an item model
}
\description{
  This function instantiates a dichotomous response model.
  The discrimination prior defaults to the lognormal
  distribution with \code{meanlog=0} and \code{sdlog=.5}.
  The guessing prior is the beta distribution. See the
  source code for details.  For discussion on the choice of
  these Bayesian priors see Baker & Kim (2004, pp.
  187-188).
}
\details{
  It is not yet possible to further customize the Bayesian
  priors. The API will change before the 1.0 release.
}
\references{
  Baker & Kim (2004). Item Response Theory: Parameter
  Estimation Techniques. Marcel Dekker, Inc.
}


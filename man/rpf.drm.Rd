\name{rpf.drm}
\alias{rpf.drm}
\title{Create a dichotomous response model and associated hyperparameters.}
\usage{
  rpf.drm(numChoices = 5, dimensions = 1, D = 1,
    multidimensional)
}
\arguments{
  \item{numChoices}{the number of choices in the question}

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


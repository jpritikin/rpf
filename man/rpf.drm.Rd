\name{rpf.drm}
\alias{rpf.drm}
\title{Create a dichotomous response model and associated hyperparameters.}
\usage{
  rpf.drm(numChoices = 5, dimensions = 1, multidimensional)
}
\arguments{
  \item{numChoices}{the number of choices in the question}

  \item{dimensions}{the number of dimensions}

  \item{multidimensional}{whether to use a multidimensional
  model. Defaults to \code{TRUE} when \code{dimensions>1}
  and \code{FALSE} when \code{dimensions==1}.}
}
\value{
  an item model
}
\description{
  Bayesian priors are only used to generate plausible
  random parameters. For discussion on the choice of these
  priors see Baker & Kim (2004, pp. 187-188).
}
\references{
  Baker & Kim (2004). Item Response Theory: Parameter
  Estimation Techniques. Marcel Dekker, Inc.
}


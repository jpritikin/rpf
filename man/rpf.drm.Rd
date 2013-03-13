\name{rpf.drm}
\alias{rpf.drm}
\title{Create a dichotomous response model and associated hyperparameters.}
\usage{
  rpf.drm(numChoices = 5, dimensions = 1, multidimensional,
    a.prior.sdlog = 0.5, poor = FALSE)
}
\arguments{
  \item{numChoices}{the number of choices in the question}

  \item{dimensions}{the number of dimensions}

  \item{multidimensional}{whether to use a multidimensional
  model. Defaults to \code{TRUE} when \code{dimensions>1}
  and \code{FALSE} when \code{dimensions==1}.}

  \item{a.prior.sdlog}{under construction}

  \item{poor}{if TRUE, use the traditional parameterization
  of the 1d model instead of the slope-intercept
  parameterization}
}
\value{
  an item model
}
\description{
  For discussion on the choice of priors see Cai, Yang, and
  Hansen (2011, p. 246).
}
\references{
  Cai, L., Yang, J. S., & Hansen, M. (2011). Generalized
  Full-Information Item Bifactor Analysis.
  \emph{Psychological Methods, 16}(3), 221-248.
}


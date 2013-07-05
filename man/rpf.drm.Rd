\name{rpf.drm}
\alias{rpf.drm}
\title{Create a dichotomous response model}
\usage{
  rpf.drm(factors = 1, multidimensional = TRUE,
    poor = FALSE)
}
\arguments{
  \item{factors}{the number of factors}

  \item{multidimensional}{whether to use a multidimensional
  model. Defaults to \code{TRUE} when \code{factors>1} and
  \code{FALSE} when \code{factors==1}.}

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


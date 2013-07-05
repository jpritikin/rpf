\name{rpf.ot2000.chisq1}
\alias{rpf.ot2000.chisq1}
\title{Compute S-Chi-squared fit statistic for 1 item}
\usage{
  rpf.ot2000.chisq1(spec, param, free, item, observed,
    quad = NULL)
}
\arguments{
  \item{spec}{a list of item specifications}

  \item{param}{item paramters}

  \item{free}{a matrix of the same shape as \code{param}
  indicating whether the parameter is free (TRUE) or fixed}

  \item{item}{the item of interest}

  \item{observed}{a matrix of observed raw scores by the
  outcome of the item of interest}

  \item{quad}{the quadrature rule (default is equally
  spaced intervals with 49 points)}
}
\description{
  Implements the Kang & Chen (2007) polytomous extension to
  S-Chi-squared statistic of Orlando & Thissen (2000).
  Fails in the presence of missing data.
}
\details{
  NOTE: IRTPRO 2.1 uses an equal interval quadrature rule
  by default.

  WARNING: The algorithm for collapsing low-count cells has
  not been tested thoroughly.
}
\references{
  Kang, T. and Chen, T. T. (2007). An investigation of the
  performance of the generalized S-Chisq item-fit index for
  polytomous IRT models. ACT Research Report Series.

  Orlando, M. and Thissen, D. (2000). Likelihood-Based
  Item-Fit Indices for Dichotomous Item Response Theory
  Models. \emph{Applied Psychological Measurement, 24}(1),
  50-64.
}


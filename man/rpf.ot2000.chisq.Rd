\name{rpf.ot2000.chisq}
\alias{rpf.ot2000.chisq}
\title{Compute S-Chi-squared fit statistic for a set of items}
\usage{
  rpf.ot2000.chisq(spec, param, free, data, quad = NULL)
}
\arguments{
  \item{spec}{a list of item specifications}

  \item{param}{item paramters in columns}

  \item{free}{a matrix of the same shape as \code{param}
  indicating whether the parameter is free (TRUE) or fixed}

  \item{data}{a data frame or matrix of response patterns,
  1 per row}

  \item{quad}{the quadrature rule (default is G-H with 49
  points)}
}
\value{
  a list of output from \code{\link{rpf.ot2000.chisq1}}
}
\description{
  Runs \code{\link{rpf.ot2000.chisq1}} for every item and
  accumulates the results.
}


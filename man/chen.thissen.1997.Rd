\name{chen.thissen.1997}
\alias{chen.thissen.1997}
\title{Computes local dependence indices for all pairs of items}
\usage{
  chen.thissen.1997(grp, data, inames = NULL, qwidth = 6,
    qpoints = 49, method = "rms")
}
\arguments{
  \item{grp}{a list with the spec, param, mean, and cov
  describing the group}

  \item{data}{data}

  \item{inames}{a subset of items to examine}

  \item{qwidth}{quadrature width}

  \item{qpoints}{number of equally spaced quadrature
  points}

  \item{method}{method to use to calculate P values. The
  default ("rms") uses the root mean square statistic (see
  \code{\link{ptw2011.gof.test}}). To obtain the
  traditional Pearson X^2 statistic, use method="pearson".}
}
\value{
  a lower triangular matrix of log P values with the sign
  determined by relative association between the observed
  and expected tables (see \code{\link{ordinal.gamma}})
}
\description{
  Item Factor Analysis makes two assumptions: (1) that the
  latent distribution is reasonably approximated by the
  multivariate Normal and (2) that items are conditionally
  independent. This test examines the second assumption.
  The presence of locally dependent items can inflate the
  precision of estimates causing the test to seem more
  accurate than it really is.
}
\details{
  Statically significant entries suggest that the item pair
  has local dependence. Since log(.01)=-4.6, an absolute
  magitude of 5 is a reasonable cut-off. Positive entries
  indicate that the two items are more correlated than
  expected. These items may share an unaccounted for latent
  dimension. Consider a redesign of the items or the use of
  testlets for scoring. Negative entries indicate that the
  two items are less correlated than expected.

  TODO: Optimize integration for the two-tier restriction.
}
\references{
  Chen, W.-H. & Thissen, D. (1997). Local dependence
  indexes for item pairs using Item Response Theory.
  \emph{Journal of Educational and Behavioral Statistics,
  22}(3), 265-289.

  Wainer, H. & Kiely, G. L. (1987). Item clusters and
  computerized adaptive testing: A case for testlets.
  \emph{Journal of Educational measurement, 24}(3),
  185--201.
}


\name{rpf.GaussHermiteData}
\alias{rpf_GaussHermiteData}
\alias{rpf.GaussHermiteData}
\title{Compute Gauss-Hermite quadrature rule}
\usage{
  rpf.GaussHermiteData(n)
}
\arguments{
  \item{n}{Order of Gauss-Hermite rule to compute (number
  of nodes)}
}
\value{
  A list containing the node positions (x) and the
  quadrature weights (w) for the requested rule.
}
\description{
  Computes Gauss-Hermite quadrature rule of requested order
  using the Golub-Welsch algorithm. This is very fast and
  numerically stable. It can handle quadrature of order
  1000+.
}
\author{
  Alexander W Blocker <ablocker@gmail.com>
}
\references{
  Golub, G. H. and Welsch, J. H. (1969). Calculation of
  Gauss Quadrature Rules. Mathematics of Computation 23
  (106): 221-230
}


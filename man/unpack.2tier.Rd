\name{unpack.2tier}
\alias{unpack.2tier}
\title{Unpack a two-tier model}
\usage{
unpack.2tier(grp)
}
\arguments{
  \item{grp}{a list containing the design matrix, spec,
  param, mean, cov, scores, and data}
}
\value{
the same model without a design matrix
}
\description{
A design matrix can be used to compress a bifactor or
two-tier model into primary factors and a partition of
items into specific factors. This function unpacks the
design matrix such that the specific factors are given
their own separate loading across all the items.
}
\examples{
# Need example here, TODO
}


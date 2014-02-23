\name{rpf.sample}
\alias{rpf.sample}
\title{Randomly sample response patterns given a list of items}
\usage{
rpf.sample(theta, items, params, design, prefix = "i", mean = NULL,
  cov = NULL)
}
\arguments{
  \item{theta}{either a vector (for 1 dimension) or a
  matrix (for >1 dimension) of person abilities or the
  number of response patterns to generate randomly}

  \item{items}{a list of item models}

  \item{params}{a list or matrix of item parameters. If
  omitted, random item parameters are generated for each
  item model.}

  \item{design}{a matrix assigning person abilities to item
  factors}

  \item{prefix}{Column names are taken from param or items.
  If no column names are available, some will be generated
  using the given prefix.}

  \item{mean}{mean vector of latent distribution
  (optional)}

  \item{cov}{covariance matrix of latent distribution
  (optional)}
}
\value{
Returns a data frame of response patterns
}
\description{
Returns a random sample of response patterns given a list
of item models and parameters.
}
\details{
The design matrix can accomodate more person abilities than
item dimension. Refer to Cai (2010) for design matrix
examples.

TODO: Add restrictions to design matrix to match
restrictions imposed by Cai (2010).
}
\examples{
# 1 dimensional items
i1 <- rpf.drm()
i1.p <- rpf.rparam(i1)
i2 <- rpf.nrm(outcomes=3)
i2.p <- rpf.rparam(i2)
rpf.sample(5, list(i1,i2), list(i1.p, i2.p))

# multidimensional items
numItems <- 4
items <- vector("list", numItems)
correct <- vector("list", numItems)

i1 <- rpf.drm(factors=2)
i2 <- rpf.drm(factors=1, multidimensional=TRUE)

for (ix in 1:(numItems-1)) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1)
}
items[[4]] <- i2
correct[[4]] <- rpf.rparam(i2)

design <- matrix(c(1, 1, 1, 1,
                   2, 2, 3, NA), nrow=2, byrow=TRUE)
rpf.sample(10, items, correct, design)
}
\references{
Cai, L. (2010). A two-tier full-information item factor
analysis model with applications. \emph{Psychometrika, 75},
581-612.
}
\seealso{
\code{\link{sample}}
}


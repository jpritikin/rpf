\name{rpf.sample}
\alias{rpf.sample}
\title{Randomly sample response patterns given a list of items}
\usage{
  rpf.sample(theta, items, params, design)
}
\arguments{
  \item{theta}{either a vector of trait abilities or the
  number of abilities to draw from N(0,1)}

  \item{items}{a list of item models}

  \item{params}{a list of item parameters. If omitted,
  random item parameters are generated for each item
  model.}

  \item{design}{assigns person abilities to item
  dimensions}
}
\value{
  Returns a matrix of response patterns
}
\description{
  Returns a random sample of response patterns given a list
  of item models and parameters.
}
\examples{
# 1 dimensional items
i1 <- rpf.drm()
i1.p <- rpf.rparam(i1)
i2 <- rpf.gpcm(numOutcomes=3)
i2.p <- rpf.rparam(i2)
rpf.sample(5, list(i1,i2), list(i1.p, i2.p))

# multidimensional items
numItems <- 4
items <- vector("list", numItems)
correct <- vector("list", numItems)

i1 <- rpf.drm(dimensions=2)
i2 <- rpf.drm(dimensions=1, multidimensional=TRUE)

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
\seealso{
  \code{\link{sample}}
}


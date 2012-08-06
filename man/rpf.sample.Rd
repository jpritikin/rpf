\name{rpf.sample}
\alias{rpf.sample}
\title{Randomly sample response patterns given a list of items}
\usage{
  rpf.sample(theta, items, params)
}
\arguments{
  \item{theta}{either a vector of trait abilities or the
  number of abilities to draw from N(0,1)}

  \item{items}{a list of item models}

  \item{params}{a list of item parameters. If omitted,
  random item parameters are generated for each item
  model.}
}
\value{
  Returns a matrix of response patterns
}
\description{
  Returns a random sample of response patterns given a list
  of item models and parameters.
}
\examples{
i1 <- rpf.drm()
i1.p <- rpf.rparam(i1)
i2 <- rpf.gpcm(numOutcomes=3)
i2.p <- rpf.rparam(i2)
rpf.sample(5, list(i1,i2), list(i1.p, i2.p))
}
\seealso{
  \code{\link{sample}}
}


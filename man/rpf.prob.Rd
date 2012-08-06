\docType{methods}
\name{rpf.prob}
\alias{rpf.prob}
\alias{rpf.prob,rpf.drm,numeric,numeric-method}
\alias{rpf.prob,rpf.gpcm,numeric,numeric-method}
\title{Map an item model, item parameters, and person trait score into a
probability vector}
\arguments{
  \item{m}{an item model}

  \item{param}{item parameters}

  \item{theta}{a single person's trait score}
}
\value{
  a vector of probabilities
}
\description{
  Map an item model, item parameters, and person trait
  score into a probability vector
}
\examples{
i1 <- rpf.drm()
i1.p <- rpf.rparam(i1)
rpf.prob(i1, c(i1.p), -1)   # low trait score
rpf.prob(i1, c(i1.p), 0)    # average trait score
rpf.prob(i1, c(i1.p), 1)    # high trait score
}
\seealso{
  See \code{\link{rpf.drm}} and \code{\link{rpf.gpcm}} to
  create item models.

  See \code{\link{rpf.rparam}} to create item parameters.
}


\docType{methods}
\name{rpf.prob}
\alias{rpf.prob}
\alias{rpf.prob,rpf.1dim.drm,numeric,numeric-method}
\alias{rpf.prob,rpf.1dim.gpcm,numeric,numeric-method}
\alias{rpf.prob,rpf.mdim.drm,numeric,matrix-method}
\alias{rpf.prob,rpf.mdim.gpcm,numeric,matrix-method}
\alias{rpf.prob,rpf.mdim.mcm,numeric,matrix-method}
\alias{rpf.prob,rpf.mdim.nrm,numeric,matrix-method}
\title{Map an item model, item parameters, and person trait score into a
probability vector}
\arguments{
  \item{m}{an item model}

  \item{param}{item parameters}

  \item{theta}{the trait score(s)}
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
rpf.prob(i1, c(i1.p), c(0,1))    # average and high trait score
}


\docType{methods}
\name{rpf.logprob}
\alias{rpf.logprob}
\alias{rpf.logprob,rpf.1dim.drm,numeric,numeric-method}
\alias{rpf.logprob,rpf.1dim.gpcm,numeric,numeric-method}
\alias{rpf.logprob,rpf.mdim.drm,numeric,matrix-method}
\title{Map an item model, item parameters, and person trait score into a
probability vector}
\arguments{
  \item{m}{an item model}

  \item{param}{item parameters}

  \item{theta}{the trait score(s)}
}
\value{
  a vector of probabilities. For dichotomous items,
  probabilities are returned in the order incorrect,
  correct. Although redundent, both incorrect and correct
  probabilities are returned in the dichotomous case for
  API consistency with polytomous item models.
}
\description{
  Map an item model, item parameters, and person trait
  score into a probability vector
}
\examples{
i1 <- rpf.gpcm()
i1.p <- rpf.rparam(i1)
rpf.logprob(i1, c(i1.p), -1)   # low trait score
rpf.logprob(i1, c(i1.p), c(0,1))    # average and high trait score
}


\docType{methods}
\name{rpf.prob}
\alias{rpf.prob}
\alias{rpf.prob,rpf.1dim.grm,numeric,numeric-method}
\alias{rpf.prob,rpf.1dim,numeric,matrix-method}
\alias{rpf.prob,rpf.1dim,numeric,numeric-method}
\alias{rpf.prob,rpf.base,data.frame,numeric-method}
\alias{rpf.prob,rpf.base,matrix,matrix-method}
\alias{rpf.prob,rpf.base,matrix,numeric-method}
\alias{rpf.prob,rpf.mdim.gpcm,numeric,matrix-method}
\alias{rpf.prob,rpf.mdim.grm,numeric,matrix-method}
\alias{rpf.prob,rpf.mdim.grm,numeric,numeric-method}
\alias{rpf.prob,rpf.mdim.mcm,numeric,matrix-method}
\alias{rpf.prob,rpf.mdim.nrm,numeric,matrix-method}
\alias{rpf.prob,rpf.mdim,numeric,matrix-method}
\alias{rpf.prob,rpf.mdim,numeric,numeric-method}
\alias{rpf_prob_wrapper}
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
  In some cases, this function is implemented in terms of
  \code{\link{rpf.logprob}}.
}
\examples{
i1 <- rpf.drm()
i1.p <- rpf.rparam(i1)
rpf.prob(i1, c(i1.p), -1)   # low trait score
rpf.prob(i1, c(i1.p), c(0,1))    # average and high trait score
}


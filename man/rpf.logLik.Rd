\docType{methods}
\name{rpf.logLik}
\alias{rpf.logLik}
\alias{rpf.logLik,rpf.1dim.drm,numeric-method}
\alias{rpf.logLik,rpf.1dim.gpcm,numeric-method}
\alias{rpf.logLik,rpf.mdim.drm,numeric-method}
\alias{rpf.logLik,rpf.mdim.gpcm,numeric-method}
\alias{rpf.logLik,rpf.mdim.mcm,numeric-method}
\alias{rpf.logLik,rpf.mdim.nrm,numeric-method}
\title{Log likelihood of item parameters with respect to Bayesian prior}
\arguments{
  \item{m}{an item model}

  \item{param}{item parameters}
}
\value{
  a log likelihood (not -2 * log likelihood)
}
\description{
  This function calculates the log likelihood of the item
  parameters with respect to the item model's Bayesian
  prior. This function is typically used by model fitting
  algorithms.
}
\examples{
i1 <- rpf.drm()
i1.p <- rpf.rparam(i1)
rpf.logLik(i1, c(i1.p))
}


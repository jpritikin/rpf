\docType{methods}
\name{rpf.rparam}
\alias{rpf.rparam}
\alias{rpf.rparam,rpf.drm-method}
\alias{rpf.rparam,rpf.gpcm-method}
\title{Generates item parameters}
\arguments{
  \item{m}{an item model}
}
\value{
  item parameters
}
\description{
  This function generates random item parameters taking the
  Bayesian priors into account.
}
\examples{
i1 <- rpf.drm()
rpf.rparam(i1)
}


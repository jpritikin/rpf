\docType{methods}
\name{rpf.rparam}
\alias{rpf.rparam}
\alias{rpf.rparam,rpf.1dim.drm-method}
\alias{rpf.rparam,rpf.1dim.graded-method}
\alias{rpf.rparam,rpf.mdim.drm-method}
\alias{rpf.rparam,rpf.mdim.graded-method}
\alias{rpf.rparam,rpf.mdim.mcm-method}
\alias{rpf.rparam,rpf.mdim.nrm-method}
\title{Generates item parameters}
\arguments{
  \item{m}{an item model}
}
\value{
  item parameters
}
\description{
  This function generates random item parameters.
}
\examples{
i1 <- rpf.drm()
rpf.rparam(i1)
}


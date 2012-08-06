\docType{methods}
\name{rpf.startingParam}
\alias{rpf.startingParam}
\alias{rpf.startingParam,rpf.drm-method}
\alias{rpf.startingParam,rpf.gpcm-method}
\title{Initial item parameters}
\arguments{
  \item{m}{an item model}
}
\value{
  item parameters
}
\description{
  This function generates item parameters suitable for use
  as initial starting values in an IRT fitting algorithm.
}
\examples{
i1 <- rpf.drm()
rpf.startingParam(i1)
}


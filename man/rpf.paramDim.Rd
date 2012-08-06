\docType{methods}
\name{rpf.paramDim}
\alias{rpf.paramDim}
\alias{rpf.paramDim,rpf.drm-method}
\alias{rpf.paramDim,rpf.gpcm-method}
\title{Dimensions of an item's parameters}
\arguments{
  \item{m}{an item model}
}
\value{
  a vector of dimension sizes
}
\description{
  This method is available in case you want to dimension a
  matrix to exactly the right size to hold item parameters.
}
\examples{
i1 <- rpf.drm()
i2 <- rpf.gpcm(numOutcomes=5)
apply(sapply(list(i1,i2), rpf.paramDim), 1, max)
}


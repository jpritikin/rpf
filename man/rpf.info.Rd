\docType{methods}
\name{rpf.info}
\alias{rpf.info}
\alias{rpf.info,rpf.1dim.drm,numeric,numeric-method}
\title{Map an item model, item parameters, and person trait score into a
information vector}
\arguments{
  \item{m}{an item model}

  \item{param}{item parameters}

  \item{theta}{the trait score(s)}
}
\value{
  Fisher information
}
\description{
  Map an item model, item parameters, and person trait
  score into a information vector
}
\examples{
i1 <- rpf.drm()
i1.p <- c(.6,1,.1)
theta <- seq(0,3,.05)
plot(theta, rpf.info(i1, i1.p, theta), type="l")
}


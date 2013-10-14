\name{rpf.dLL}
\alias{rpf.dLL}
\alias{rpf.dLL,rpf.base,numeric,numeric,numeric-method}
\alias{rpf_dLL_wrapper}
\title{Item parameter derivatives}
\arguments{
  \item{m}{item model}

  \item{param}{item parameters}

  \item{where}{location in the latent space}

  \item{weight}{per outcome weights (typically derived by
  observation)}
}
\value{
  first and second order partial derivatives of the log
  likelihood evaluated at \code{where}. For p parameters,
  the first p values are the first derivative and the next
  p(p+1)/2 columns are the lower triangle of the second
  derivative.
}
\description{
  Evaluate the partial derivatives of the log likelihood
  with respect to each parameter at \code{where} with
  \code{weight}.
}
\seealso{
  The numDeriv package.
}


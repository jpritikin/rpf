\name{rpf.gradient}
\alias{rpf.gradient}
\alias{rpf.gradient,rpf.base,numeric,numeric,numeric-method}
\alias{rpf_gradient_wrapper}
\title{Item parameter gradients}
\arguments{
  \item{m}{item model}

  \item{param}{item parameters}

  \item{where}{location in the latent space}

  \item{weight}{per outcome weights (typically derived by
  observation)}
}
\value{
  derivative of the log likelihood with respect to each
  parameter evaluated at \code{where}
}
\description{
  Evaluate the partial derivatives of the log likelihood
  with respect to each parameter at \code{where} with
  \code{weight}.
}


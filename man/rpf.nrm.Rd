\name{rpf.nrm}
\alias{rpf.nrm}
\title{Create a nominal response model}
\usage{
  rpf.nrm(outcomes = 3, factors = 1, T.a = "trend",
    T.c = "trend")
}
\arguments{
  \item{outcomes}{The number of choices available}

  \item{factors}{the number of factors}

  \item{T.a}{the T matrix for slope parameters}

  \item{T.c}{the T matrix for intercept parameters}
}
\value{
  an item model
}
\description{
  This function instantiates a nominal response model. The
  T matrix must be an invertible square matrix of dimension
  outcomes-1. As a shortcut, either T matrix can be
  specified as "trend" for a Fourier basis or as "id" for
  an identity basis.
}
\examples{
# typical parameterization for the Generalized Partial Credit Model
gpcm <- function(outcomes) rpf.nrm(outcomes, T.c=lower.tri(diag(outcomes-1),TRUE) * -1)
}
\references{
  Thissen, D., Cai, L., & Bock, R. D. (2010). The Nominal
  Categories Item Response Model. In M. L. Nering & R.
  Ostini (Eds.), \emph{Handbook of Polytomous Item Response
  Theory Models} (pp. 43--75). Routledge.
}


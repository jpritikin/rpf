\name{rpf.nrm}
\alias{rpf.nrm}
\title{Create a nominal response model}
\usage{
rpf.nrm(outcomes = 3, factors = 1, T.a = "trend", T.c = "trend")
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
This function instantiates a nominal response model.
}
\details{
The transformation matrices T.a and T.c are chosen by the
analyst and not estimated.  The T matrices must be
invertible square matrices of dimension outcomes-1. As a
shortcut, either T matrix can be specified as "trend" for a
Fourier basis or as "id" for an identity basis. The
response probability function is

\deqn{a_k = T_a \alpha_k} \deqn{c_k = T_c \gamma_k}
\deqn{\mathrm P(\mathrm{pick}=k|a,a_k,c_k) = C\
\frac{1}{1+\exp(-(a \theta a_k + c_k))}}

where \eqn{a_k} and \eqn{c_k} are the result of multiplying
two vectors of free parameters \eqn{\alpha} and
\eqn{\gamma} by fixed matrices \eqn{T_a} and \eqn{T_c},
respectively; \eqn{a_0} and \eqn{c_0} are fixed to 0 for
identification; and \eqn{C} is a normalizing constant to
ensure that \eqn{\sum_k \mathrm P(\mathrm{pick}=k) = 1}.
}
\examples{
# typical parameterization for the Generalized Partial Credit Model
gpcm <- function(outcomes) rpf.nrm(outcomes, T.c=lower.tri(diag(outcomes-1),TRUE) * -1)
}
\references{
Thissen, D., Cai, L., & Bock, R. D. (2010). The Nominal
Categories Item Response Model. In M. L. Nering & R. Ostini
(Eds.), \emph{Handbook of Polytomous Item Response Theory
Models} (pp. 43--75). Routledge.
}


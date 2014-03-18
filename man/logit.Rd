\name{logit}
\alias{logit}
\title{Transform from [0,1] to the reals}
\usage{
logit(prob)
}
\arguments{
  \item{prob}{a number between 0 and 1}
}
\description{
The logit function is a standard transformation from [0,1]
(such as a probability) to the real number line.
}
\examples{
logit(.5)  # 0
logit(.25) # -1.098
logit(0)   # -Inf
}


\name{write.flexmirt}
\alias{write.flexmirt}
\title{Write a flexMIRT PRM file}
\usage{
write.flexmirt(groups, file = NULL, fileEncoding = "")
}
\arguments{
  \item{groups}{a list of groups each with items and latent
  parameters}

  \item{file}{the destination file name}

  \item{fileEncoding}{how to encode the text file
  (optional)}
}
\description{
Formats item parameters in the way that flexMIRT expects to
read them. Use \code{\link{read.flexmirt}} to see what
shape the groups parameter of this function should take.
}
\details{
NOTE: Support for the graded response model is not
complete.
}


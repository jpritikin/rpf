\name{rpf.grm}
\alias{rpf.grm}
\title{Create a graded response model}
\usage{
  rpf.grm(outcomes = 2, factors = 1,
    multidimensional = TRUE)
}
\arguments{
  \item{outcomes}{The number of choices available}

  \item{factors}{the number of factors}

  \item{multidimensional}{whether to use a multidimensional
  model. Defaults to \code{TRUE}.}
}
\value{
  an item model
}
\description{
  This function instantiates a graded response model.
}
\details{
  The graded response model was designed for a item with a
  series of dependent parts where a higher score implies
  that easier parts of the item were surmounted. If there
  is any chance your polytomous item has independent parts
  then consider \code{\link{rpf.nrm}}. If your categories
  cannot cross then the graded response model provides a
  little more information than the nominal model. Stronger
  a priori assumptions offer provide more power at the cost
  of flexibility.
}


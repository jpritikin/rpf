\name{rpf.1dim.fit}
\alias{rpf.1dim.fit}
\title{Calculate item and person fit statistics}
\usage{
  rpf.1dim.fit(spec, params, responses, scores, margin,
    na.rm = TRUE, wh.exact = TRUE)
}
\arguments{
  \item{spec}{list of item models}

  \item{params}{data frame of item parameters, 1 per row}

  \item{responses}{persons in rows and items in columns}

  \item{scores}{model derived person scores}

  \item{margin}{for people 1, for items 2}

  \item{na.rm}{remove NAs (default TRUE)}

  \item{wh.exact}{whether to use the exact Wilson-Hilferty
  transformation (default TRUE)}
}
\description{
  Exact distributional properties of these statistics are
  unknown (Masters & Wright, 1997, p. 112).  For details on
  the calculation, refer to Wright & Masters (1982, p.
  100).
}
\details{
  The Wilson-Hilferty transformation is biased for less
  than 25 items. To adjust Z scores for fewer items use
  wh.exact=FALSE.
}
\references{
  Masters, G. N. & Wright, B. D. (1997). The Partial Credit
  Model. In W. van der Linden & R. K. Kambleton (Eds.),
  \emph{Handbook of modern item response theory} (pp.
  101-121). Springer.

  Wilson, E. B., & Hilferty, M. M. (1931). The distribution
  of chi-square. \emph{Proceedings of the National Academy
  of Sciences of the United States of America,} 17,
  684-688.

  Wright, B. D. & Masters, G. N. (1982). \emph{Rating Scale
  Analysis.} Chicago: Mesa Press.
}


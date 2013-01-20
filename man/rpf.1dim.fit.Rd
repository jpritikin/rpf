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
  For details on the calculation, see Embretson & Reise
  (2000, pp. 237-238) or Wright & Masters (1982, p. 100).
}
\details{
  The Wilson-Hilferty transformation is biased for less
  than 25 items. To adjust Z scores for fewer items use
  wh.exact=FALSE.
}
\references{
  Embretson, S. E. & Reise S. P. (2000) \emph{Item response
  theory for psychologists.} Lawrence Erlbaum.

  Wilson, E. B., & Hilferty, M. M. (1931). The distribution
  of chi-square. \emph{Proceedings of the National Academy
  of Sciences of the United States of America,} 17,
  684-688.

  Wright, B. D. & Masters, G. N. (1982). \emph{Rating Scale
  Analysis.} Chicago: Mesa Press.
}


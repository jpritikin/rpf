##' Create a generalized partial credit model and associated hyperparameters.
##'
##' @param outcomes The number of choices available
##' @param factors the number of factors
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{TRUE} when \code{factors>1} and
##' \code{FALSE} when \code{factors==1}.
##' @return an item model
##' @export
##' @references Baker & Kim (2004). \emph{Item Response Theory: Parameter
##' Estimation Techniques.} Marcel Dekker, Inc.
##'
##' Masters, G. N. (1982). A rasch model for partial credit scoring.
##' \emph{Psychometrika, 47}(2), 149-174.
##'
##' Muraki, Eiji. (1992). A generalized partial credit model:
##' Application of an EM algorithm. \emph{Applied Psychological Measurement,
##' 16}, 159-176. doi:10.1177/014662169201600206
rpf.gpcm <- function(outcomes=2, factors=1, multidimensional) {
  if (missing(multidimensional)) {
    multidimensional <- factors > 1
  }
  if (!multidimensional && factors > 1) {
    stop("More than 1 dimension must use a multidimensional model")
  }
  if (!multidimensional) {
    m <- new("rpf.1dim.gpcm",
        outcomes=outcomes,
        factors=1)
    m@spec <- c(rpf.id_of('gpcm1'), m@outcomes, m@factors)
    m
  } else {
    m <- new("rpf.mdim.gpcm",
        outcomes=outcomes,
        factors=factors)
#    m@spec <- c(rpf.id_of('gpcm'), m@outcomes, m@factors)
    m
  }
}

### 1dim

### mdim

setMethod("rpf.prob", signature(m="rpf.mdim.gpcm", param="numeric",
                                theta="matrix"),
          function(m, param, theta) {
            a <- param[1:m@factors]
            b <- param[-1:-m@factors]
            tri <- lower.tri(matrix(NA, m@outcomes, m@outcomes))
            numPersons <- dim(theta)[1]
            out <- array(dim=c(numPersons, m@outcomes))
            for (px in 1:numPersons) {
                k <- exp(apply(c(0, -((theta[px,] %*% a) + b)) * tri, c(2),sum))
                out[px,] <- k / sum(k)
            }
            return(out)
          })

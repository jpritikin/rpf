##' Create a generalized partial credit model and associated hyperparameters.
##'
##' @param numOutcomes The number of choices available
##' @param dimensions the number of dimensions
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{TRUE} when \code{dimensions>1} and
##' \code{FALSE} when \code{dimensions==1}.
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
rpf.gpcm <- function(numOutcomes=2, dimensions=1, multidimensional) {
  if (missing(multidimensional)) {
    multidimensional <- dimensions > 1
  }
  if (!multidimensional && dimensions > 1) {
    stop("More than 1 dimension must use a multidimensional model")
  }
  if (!multidimensional) {
    m <- new("rpf.1dim.gpcm",
        numOutcomes=numOutcomes,
        dimensions=1,
        numParam=numOutcomes,
        a.prior.sdlog=.5)
    m@spec <- c(rpf.id_of('gpcm1'), m@numOutcomes, m@dimensions, m@a.prior.sdlog)
    m
  } else {
    new("rpf.mdim.gpcm",
        numOutcomes=numOutcomes,
        dimensions=dimensions,
        numParam=dimensions + numOutcomes - 1,
        a.prior.sdlog=.5)
  }
}

### 1dim

### mdim

setMethod("rpf.prob", signature(m="rpf.mdim.gpcm", param="numeric",
                                theta="matrix"),
          function(m, param, theta) {
            a <- param[1:m@dimensions]
            b <- param[-1:-m@dimensions]
            tri <- lower.tri(matrix(NA, m@numOutcomes, m@numOutcomes))
            numPersons <- dim(theta)[1]
            out <- array(dim=c(numPersons, m@numOutcomes))
            for (px in 1:numPersons) {
                k <- exp(apply(c(0, -((theta[px,] %*% a) + b)) * tri, c(2),sum))
                out[px,] <- k / sum(k)
            }
            return(out)
          })

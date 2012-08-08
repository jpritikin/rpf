##' Create a generalized partial credit model and associated hyperparameters.
##'
##' This function instantiates a generalized partial credit model. The
##' discrimination prior defaults to the lognormal distribution with
##' \code{meanlog=0} and \code{sdlog=.5}.
##'
##' It is not yet possible to further customize the Bayesian
##' priors. The API will change before the 1.0 release.
##' 
##' @param numOutcomes The number of choices available
##' @param dimensions the number of dimensions
##' @param D defaults to 1 or pass in \code{\link{rpf.ogive}}
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{TRUE} when \code{dimensions>1} and
##' \code{FALSE} when \code{dimensions==1}.
##' @return an item model
##' @export
##' @references Baker & Kim (2004). Item Response Theory: Parameter
##' Estimation Techniques. Marcel Dekker, Inc.
rpf.gpcm <- function(numOutcomes=2, dimensions=1, D=1, multidimensional) {
  if (missing(multidimensional)) {
    multidimensional <- dimensions > 1
  }
  if (!multidimensional && dimensions > 1) {
    stop("More than 1 dimension must use a multidimensional model")
  }
  if (!multidimensional) {
    new("rpf.1dim.gpcm", numOutcomes=numOutcomes, D=D,
        dimensions=1,
        numParam=numOutcomes,
        a.prior.meanlog=0,
        a.prior.sdlog=.5)
  } else {
    new("rpf.mdim.gpcm", numOutcomes=numOutcomes, D=D,
        dimensions=dimensions,
        numParam=dimensions + numOutcomes - 1,
        a.prior.meanlog=0,
        a.prior.sdlog=.5)
  }
}

### 1dim

setMethod("rpf.prob", signature(m="rpf.1dim.gpcm", param="numeric",
                                theta="numeric"),
          function(m, param, theta) {
            a <- param[1]
            b <- param[-1]
            tri <- lower.tri(matrix(NA, m@numOutcomes, m@numOutcomes))
            out <- array(dim=c(length(theta), m@numOutcomes))
            for (px in 1:length(theta)) {
                k <- exp(apply(c(0, m@D * -a * (theta[px] - b)) * tri, c(2), sum))
                out[px,] <- k / sum(k)
            }
            return(out)
          })

setMethod("rpf.logLik", signature(m="rpf.1dim.gpcm", param="numeric"),
          function(m, param) {
            a <- param[1]
            dlnorm(a, meanlog=m@a.prior.meanlog,
                        sdlog=m@a.prior.sdlog, log=TRUE)
          })

setMethod("rpf.rparam", signature(m="rpf.1dim.gpcm"),
          function(m) {
              a <- rlnorm(1, meanlog=m@a.prior.meanlog,
                          sdlog=m@a.prior.sdlog)
              b <- sort(rnorm(m@numOutcomes-1))
              c(a=a,b=b)
          })

setMethod("rpf.startingParam", signature(m="rpf.1dim.gpcm"),
          function(m) {
            c(a=1, b=rep(0, (m@numOutcomes-1)))
          })

setMethod("rpf.getLocation", signature(m="rpf.1dim.gpcm", param="numeric"),
          function(m, param) {
            param[2:m@numOutcomes]
          })

setMethod("rpf.setLocation", signature(m="rpf.1dim.gpcm", param="numeric", loc="numeric"),
          function(m, param, loc) {
              param[2:m@numOutcomes] <- loc
              param
          })

### mdim

setMethod("rpf.prob", signature(m="rpf.mdim.gpcm", param="numeric",
                                theta="matrix"),
          function(m, param, theta) {
            a <- param[1:m@dimensions] * m@D
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

setMethod("rpf.logLik", signature(m="rpf.mdim.gpcm", param="numeric"),
          function(m, param) {
            a <- param[1:m@dimensions]
            sum(dlnorm(a, meanlog=m@a.prior.meanlog,
                       sdlog=m@a.prior.sdlog, log=TRUE))
          })

setMethod("rpf.rparam", signature(m="rpf.mdim.gpcm"),
          function(m) {
              a <- rlnorm(m@dimensions, meanlog=m@a.prior.meanlog,
                          sdlog=m@a.prior.sdlog)
              b <- sort(rnorm(m@numOutcomes-1))
              c(a=a,b=b)
          })

setMethod("rpf.startingParam", signature(m="rpf.mdim.gpcm"),
          function(m) {
            c(a=rep(1,m@dimensions), b=rep(0, (m@numOutcomes-1)))
          })

setMethod("rpf.getLocation", signature(m="rpf.mdim.gpcm", param="numeric"),
          function(m, param) {
            param[-1:-m@dimensions]
          })

setMethod("rpf.setLocation", signature(m="rpf.mdim.gpcm", param="numeric", loc="numeric"),
          function(m, param, loc) {
            param[-1:-m@dimensions] <- loc
            param
          })

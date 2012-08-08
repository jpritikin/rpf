##' Create a nominal response model and associated hyperparameters.
##'
##' This function instantiates a nominal response model. The
##' discrimination prior defaults to the lognormal distribution with
##' \code{meanlog=0} and \code{sdlog=.5}.
##' 
##' It is not yet possible to further customize the Bayesian
##' priors. The API will change before the 1.0 release.
##' 
##' @param numOutcomes The number of choices available
##' @param dimensions the number of dimensions
##' @return an item model
##' @export
##' @author Jonathan Weeks <weeksjp@@gmail.com>
rpf.nrm <- function(numOutcomes=2, dimensions=1) {
    new("rpf.mdim.nrm",
        numOutcomes=numOutcomes, dimensions=dimensions,
        numParam=(1 + dimensions) * numOutcomes,
        a.prior.meanlog=0,
        a.prior.sdlog=.5)
}

### mdim

setMethod("rpf.prob", signature(m="rpf.mdim.nrm", param="numeric",
                                theta="matrix"),
          function(m, param, theta) {
            ##   Object for the denominator in the final MMCM equation
            den <- NULL
            
            numA <- m@numOutcomes * m@dimensions
            a1 <- param[1:numA]
            b1 <- param[-1:-numA]
            
            ##   Compute the denominator
            for (k in 1:m@numOutcomes) {
              tmp <- (k-1)*m@dimensions
              tmp1 <- tmp+m@dimensions
              d <- exp((theta %*% a1[(tmp+1):tmp1])+b1[k])
              den <- cbind(den, d)
            }
            den <- apply(den,1,sum)
            
            numPersons <- dim(theta)[1]

            p <- array(dim=c(numPersons, m@numOutcomes))

            for (k in 1:m@numOutcomes) {
              tmp <- (k-1)*m@dimensions
              tmp1 <- tmp+m@dimensions
              cp <- exp((theta %*% a1[(tmp+1):tmp1])+b1[k])/den
              p[,k] <- cp
            }
            return(p)
          })

setMethod("rpf.logLik", signature(m="rpf.mdim.nrm", param="numeric"),
          function(m, param) {
            a <- param[1:m@numOutcomes * m@dimensions]
            sum(dlnorm(a, meanlog=m@a.prior.meanlog,
                       sdlog=m@a.prior.sdlog, log=TRUE))
          })

setMethod("rpf.rparam", signature(m="rpf.mdim.nrm"),
          function(m) {
              a <- rlnorm(m@numOutcomes * m@dimensions,
                          meanlog=m@a.prior.meanlog,
                          sdlog=m@a.prior.sdlog)
              b <- sort(rnorm(m@numOutcomes))
              c(a=a,b=b)
          })

setMethod("rpf.startingParam", signature(m="rpf.mdim.nrm"),
          function(m) {
            c(a=rep(1,m@numOutcomes * m@dimensions),
              b=rep(0, m@numOutcomes))
          })

setMethod("rpf.getLocation", signature(m="rpf.mdim.nrm", param="numeric"),
          function(m, param) {
            param[-1:-m@numOutcomes * m@dimensions]
          })

setMethod("rpf.setLocation", signature(m="rpf.mdim.nrm", param="numeric", loc="numeric"),
          function(m, param, loc) {
            param[-1:-m@numOutcomes * m@dimensions] <- loc
            param
          })

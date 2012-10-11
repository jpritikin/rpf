##' Create a multiple-choice response model and associated hyperparameters.
##'
##' This function instantiates a multiple-choice response model. The
##' discrimination prior defaults to the lognormal distribution with
##' \code{meanlog=0} and \code{sdlog=.5}. The guessing prior is the
##' beta distribution. See the source code for details.
##' 
##' It is not yet possible to further customize the Bayesian
##' priors. The API will change before the 1.0 release.
##' 
##' @param numOutcomes the number of possible outcomes
##' @param numChoices the number of choices available
##' @param dimensions the number of dimensions
##' @return an item model
##' @export
##' @author Jonathan Weeks <weeksjp@@gmail.com>
rpf.mcm <- function(numOutcomes=2, numChoices=5, dimensions=1) {
  guess.weight <- 20
  guessing <- (1/numChoices)
  new("rpf.mdim.mcm",
      numOutcomes=numOutcomes, dimensions=dimensions,
      numParam=numOutcomes * (dimensions + 2) - 1,
      a.prior.meanlog=0,
      a.prior.sdlog=.5,
      c.prior.alpha=guess.weight*guessing+1,
      c.prior.beta=guess.weight*(1-guessing)+1)
}

### mdim

setMethod("rpf.prob", signature(m="rpf.mdim.mcm", param="numeric",
                                theta="matrix"),
          function(m, param, theta) {
            den <- NULL
            
            a1 <- param[1:(m@dimensions*m@numOutcomes)]
            b1 <- param[(m@dimensions*m@numOutcomes+1):
                        ((m@dimensions+1)*m@numOutcomes)]
            c1 <- param[-1:-((m@dimensions+1)*m@numOutcomes)]
            
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
              if (k==1) {
                cp <- (exp((theta %*% a1[(tmp+1):tmp1])+b1[k]))/den
              } else {
                cp <- (exp((theta %*% a1[(tmp+1):tmp1])+b1[k])+c1[k-1]*(exp((theta %*% a1[1:m@dimensions])+b1[1])))/den
              }
              p[,k] <- cp
            }
            return(p)
          })

setMethod("rpf.rparam", signature(m="rpf.mdim.mcm"),
          function(m) {
              a <- rlnorm(m@numOutcomes * m@dimensions,
                          meanlog=m@a.prior.meanlog,
                          sdlog=m@a.prior.sdlog)
              b <- sort(rnorm(m@numOutcomes))
              c <- rbeta(m@numOutcomes-1, shape1=m@c.prior.alpha-2,
                         shape2=m@c.prior.beta-2)
              c(a=a,b=b,c=c)
          })

setMethod("rpf.startingParam", signature(m="rpf.mdim.mcm"),
          function(m) {
            c(a=rep(1,m@numOutcomes * m@dimensions),
              b=rep(0, m@numOutcomes),
              c=rep(0, m@numOutcomes-1))
          })

setMethod("rpf.getLocation", signature(m="rpf.mdim.mcm", param="numeric"),
          function(m, param) {
            param[(m@dimensions*m@numOutcomes+1):
                  ((m@dimensions+1)*m@numOutcomes)]
          })

setMethod("rpf.setLocation", signature(m="rpf.mdim.mcm", param="numeric", loc="numeric"),
          function(m, param, loc) {
            param[(m@dimensions*m@numOutcomes+1):
                  ((m@dimensions+1)*m@numOutcomes)] <- loc
            param
          })

##' Create a dichotomous response model and associated hyperparameters.
##'
##' This function instantiates a dichotomous response model. The
##' discrimination prior defaults to the lognormal distribution with
##' \code{meanlog=0} and \code{sdlog=.5}. The guessing prior is the
##' beta distribution. See the source code for details.  For
##' discussion on the choice of these Bayesian priors see Baker & Kim
##' (2004, pp. 187-188).
##'
##' It is not yet possible to further customize the Bayesian
##' priors. The API will change before the 1.0 release.
##' 
##' @param numChoices the number of choices in the question
##' @param dimensions the number of dimensions
##' @param D defaults to 1 or pass in \code{\link{rpf.ogive}}
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{TRUE} when \code{dimensions>1} and
##' \code{FALSE} when \code{dimensions==1}.
##' @return an item model
##' @export
##' @references Baker & Kim (2004). Item Response Theory: Parameter
##' Estimation Techniques. Marcel Dekker, Inc.
rpf.drm <- function(numChoices=5, dimensions=1, D=1, multidimensional) {
  if (missing(multidimensional)) {
    multidimensional <- dimensions > 1
  }
  if (!multidimensional && dimensions > 1) {
    stop("More than 1 dimension must use a multidimensional model")
  }
  guess.weight <- 20
  guessing <- (1/numChoices)
  if (!multidimensional) {
    new("rpf.1dim.drm", numOutcomes=2, D=D, numParam=3,
        dimensions=1,
        guessing=guessing,
        a.prior.meanlog=0,
        a.prior.sdlog=.5,
        c.prior.alpha=guess.weight*guessing+1,
        c.prior.beta=guess.weight*(1-guessing)+1)
  } else {
    new("rpf.mdim.drm", numOutcomes=2, D=D, dimensions=dimensions,
        numParam=2+dimensions,
        guessing=guessing,
        a.prior.meanlog=0,
        a.prior.sdlog=.5,
        c.prior.alpha=guess.weight*guessing+1,
        c.prior.beta=guess.weight*(1-guessing)+1)
  }
}

### 1dim

setMethod("rpf.prob", signature(m="rpf.1dim.drm", param="numeric",
                                theta="numeric"),
          function(m, param, theta) {
            a <- param[1] * m@D
            b <- param[2]
            c <- param[3]
            p <- c + (1-c)/(1+exp(-a*(theta-b)))
            cbind(1-p,p)
          })

##' @references
##' Embretson & Reise (2000, p. 184)
setMethod("rpf.info", signature(m="rpf.1dim.drm", param="numeric",
                                theta="numeric"),
          function(m, param, theta) {
            p <- rpf.prob(m, param, theta)
            a <- param[1] * m@D
            c <- param[3]
            a^2 * (p[,1]/p[,2]) * ((p[,2]-c)^2/(1-c)^2)
          })

setMethod("rpf.rparam", signature(m="rpf.1dim.drm"),
          function(m) {
            n <- 1
            c(a=rlnorm(n, meanlog=m@a.prior.meanlog,
                sdlog=m@a.prior.sdlog),
              b=rnorm(n),
              c=rbeta(n, shape1=m@c.prior.alpha-2,
                shape2=m@c.prior.beta-2))
          })

setMethod("rpf.startingParam", signature(m="rpf.1dim.drm"),
          function(m) {
            c(a=1, b=0, c=m@guessing)
          })

setMethod("rpf.getLocation", signature(m="rpf.1dim.drm", param="numeric"),
          function(m, param) {
              param[2]
          })

setMethod("rpf.setLocation", signature(m="rpf.1dim.drm", param="numeric", loc="numeric"),
          function(m, param, loc) {
              param[2] <- loc
              param
          })

### mdim

##' @author Jonathan Weeks <weeksjp@@gmail.com>
setMethod("rpf.prob", signature(m="rpf.mdim.drm", param="numeric",
                                theta="matrix"),
          function(m, param, theta) {
            a <- param[1:m@dimensions] * m@D
            p.rest <- param[-1:-m@dimensions]
            b <- p.rest[1]
            c <- p.rest[2]
            p <- c + (1-c)/(1+exp(-(theta %*% a + b)))
            cbind(1-p,p)
          })

setMethod("rpf.rparam", signature(m="rpf.mdim.drm"),
          function(m) {
            c(a=rlnorm(m@dimensions, meanlog=m@a.prior.meanlog,
                sdlog=m@a.prior.sdlog),
              b=rnorm(1),
              c=rbeta(1, shape1=m@c.prior.alpha-2,
                shape2=m@c.prior.beta-2))
          })

setMethod("rpf.startingParam", signature(m="rpf.mdim.drm"),
          function(m) {
            c(a=rep(1,m@dimensions), b=0, c=m@guessing)
          })

setMethod("rpf.getLocation", signature(m="rpf.mdim.drm", param="numeric"),
          function(m, param) {
            param[m@dimensions+1]
          })

setMethod("rpf.setLocation", signature(m="rpf.mdim.drm", param="numeric",
                                       loc="numeric"),
          function(m, param, loc) {
            param[m@dimensions+1] <- loc
            param
          })

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
##' @param D defaults to 1 or pass in the \code{\link{rpf.ogive}}
##' @param numChoices the number of alternatives in the question
##' @return an item model
##' @export
##' @references Baker & Kim (2004). Item Response Theory: Parameter
##' Estimation Techniques. Marcel Dekker, Inc.
rpf.drm <- function(D=1, numChoices=5) {
  guess.weight <- 20
  guessing <- (1/numChoices)
  new("rpf.drm", numOutcomes=2, D=D,
      guessing=guessing,
      a.prior.meanlog=0,
      a.prior.sdlog=.5,
      c.prior.alpha=guess.weight*guessing+1,
      c.prior.beta=guess.weight*(1-guessing)+1)
}

setMethod("rpf.prob", signature(m="rpf.drm", param="numeric",
                                theta="numeric"),
          function(m, param, theta) {
            a <- param[1]
            b <- param[2]
            c <- param[3]
            p <- c + (1-c)/(1+exp(-m@D*a*(theta-b)))
            cbind(1-p,p)
          })

setMethod("rpf.logLik", signature(m="rpf.drm", param="numeric"),
          function(m, param) {
            a <- param[1]
            c <- param[3]
            sum(dlnorm(a, meanlog=m@a.prior.meanlog,
                               sdlog=m@a.prior.sdlog, log=TRUE),
                     dbeta(c, shape1=m@c.prior.alpha-2,
                           shape2=m@c.prior.beta-2, log=TRUE))
          })

setMethod("rpf.paramDim", signature(m="rpf.drm"), function(m) c(1,3))

setMethod("rpf.rparam", signature(m="rpf.drm"),
          function(m) {
              n <- 1
              cbind(a=rlnorm(n, meanlog=m@a.prior.meanlog,
                       sdlog=m@a.prior.sdlog),
                b=rnorm(n),
                c=rbeta(n, shape1=m@c.prior.alpha-2,
                      shape2=m@c.prior.beta-2))
          })

setMethod("rpf.startingParam", signature(m="rpf.drm"),
          function(m) {
              cbind(a=1, b=0, c=m@guessing)
          })

setMethod("rpf.getLocation", signature(m="rpf.drm", param="numeric"),
          function(m, param) {
              param[2]
          })

setMethod("rpf.setLocation", signature(m="rpf.drm", param="numeric", loc="numeric"),
          function(m, param, loc) {
              param[2] <- loc
              param
          })

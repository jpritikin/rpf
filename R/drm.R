# For discussion of Bayesian priors, see Baker & Kim (2004, pp. 187-188)

rpf.drm <- function(D=1, numAlternatives=5) {
  guess.weight <- 20
  guessing <- (1/numAlternatives)
  new("rpf.drm", numOutcomes=2, D=D,
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
            -2 * sum(dlnorm(a, meanlog=m@a.prior.meanlog,
                               sdlog=m@a.prior.sdlog, log=TRUE),
                     dbeta(c, shape1=m@c.prior.alpha-2,
                           shape2=m@c.prior.beta-2, log=TRUE))
          })

setMethod("rpf.rparam", signature(m="rpf.drm"),
          function(m) {
              n <- 1
              cbind(a=rlnorm(n, meanlog=m@a.prior.meanlog,
                       sdlog=m@a.prior.sdlog),
                b=rnorm(n),
                c=rbeta(n, shape1=m@c.prior.alpha-2,
                      shape2=m@c.prior.beta-2))
          })

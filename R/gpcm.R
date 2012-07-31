rpf.gpcm <- function(numChoices=2, D=1) {
  new("rpf.gpcm", numChoices=numChoices, D=D,
      a.prior.meanlog=0,
      a.prior.sdlog=.5)
}

setMethod("rpf.prob", signature(m="rpf.gpcm", param="numeric",
                                theta="numeric"),
          function(m, param, theta) {
            a <- param[1]
            b <- param[-1]
            tri <- lower.tri(matrix(NA, m@numChoices, m@numChoices))
            k <- exp(apply(c(0, m@D * -a * (theta - b)) * tri, c(2), sum))
            k / sum(k)
          })

setMethod("rpf.logLik", signature(m="rpf.gpcm", param="numeric"),
          function(m, param) {
            a <- param[1]
            -2 * dlnorm(p.a, meanlog=m@a.prior.meanlog,
                        sdlog=m@a.prior.sdlog, log=TRUE)
          })

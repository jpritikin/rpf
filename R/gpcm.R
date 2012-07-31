rpf.gpcm <- function(numOutcomes=2, D=1) {
  new("rpf.gpcm", numOutcomes=numOutcomes, D=D,
      a.prior.meanlog=0,
      a.prior.sdlog=.5)
}

setMethod("rpf.prob", signature(m="rpf.gpcm", param="numeric",
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

setMethod("rpf.logLik", signature(m="rpf.gpcm", param="numeric"),
          function(m, param) {
            a <- param[1]
            -2 * dlnorm(p.a, meanlog=m@a.prior.meanlog,
                        sdlog=m@a.prior.sdlog, log=TRUE)
          })

setMethod("rpf.rparam", signature(m="rpf.gpcm"),
          function(m) {
              a <- rlnorm(1, meanlog=m@a.prior.meanlog,
                          sdlog=m@a.prior.sdlog)
              c(a=a, b=sort(rnorm(m@numOutcomes-1)))
          })

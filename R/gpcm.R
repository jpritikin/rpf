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
            dlnorm(a, meanlog=m@a.prior.meanlog,
                        sdlog=m@a.prior.sdlog, log=TRUE)
          })

setMethod("rpf.paramDim", signature(m="rpf.gpcm"), function(m) {
    c(1, m@numOutcomes)
})

setMethod("rpf.rparam", signature(m="rpf.gpcm"),
          function(m) {
              n <- 1
              a <- rlnorm(n, meanlog=m@a.prior.meanlog,
                          sdlog=m@a.prior.sdlog)
              b <- array(dim=c(n, m@numOutcomes-1),
                         data=rnorm(n * (m@numOutcomes-1)))
              if (m@numOutcomes > 2) {
                  b <- t(apply(b, 1, sort))
              }
              out <- cbind(a=a, b)
              dimnames(out)[[2]][-1] <- paste(sep='','b',1:(m@numOutcomes-1))
              return(out)
          })

setMethod("rpf.startingParam", signature(m="rpf.gpcm"),
          function(m) {
              t(as.matrix(c(a=1, b=rep(0, (m@numOutcomes-1)))))
          })

setMethod("rpf.getLocation", signature(m="rpf.gpcm", param="numeric"),
          function(m, param) {
              t(param[2:m@numOutcomes])
          })

setMethod("rpf.setLocation", signature(m="rpf.gpcm", param="numeric", loc="numeric"),
          function(m, param, loc) {
              param[2:m@numOutcomes] <- loc
              param
          })

### 1dim

setMethod("rpf.info", signature(m="rpf.1dim.graded", param="numeric", theta="numeric"),
          function(m, param, theta) {
            pr <- rpf.prob(m, param, theta)
            tbar <- apply(pr, 1, function(r) sum(r * 1:m@numOutcomes))
            ret <- numeric(length(theta))
            for (rx in 1:length(theta)) {
              ret[rx] <- sum((1:m@numOutcomes - tbar[rx])^2 * pr[rx,])
            }
            ret * param[1]^2
          })

setMethod("rpf.rparam", signature(m="rpf.1dim.graded"),
          function(m) {
              a <- rlnorm(1, meanlog=0, sdlog=m@a.prior.sdlog)
              b <- sort(rnorm(m@numOutcomes-1))
              c(a=a,b=b)
          })

### mdim

setMethod("rpf.rparam", signature(m="rpf.mdim.graded"),
          function(m) {
              a <- rlnorm(m@dimensions, meanlog=0, sdlog=m@a.prior.sdlog)
              b <- sort(rnorm(m@numOutcomes-1))
              c(a=a,b=b)
          })

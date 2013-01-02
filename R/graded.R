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
              a <- rlnorm(1, meanlog=m@a.prior.meanlog,
                          sdlog=m@a.prior.sdlog)
              b <- sort(rnorm(m@numOutcomes-1))
              c(a=a,b=b)
          })

setMethod("rpf.startingParam", signature(m="rpf.1dim.graded"),
          function(m) {
            c(a=1, b=rep(0, (m@numOutcomes-1)))
          })

setMethod("rpf.getLocation", signature(m="rpf.1dim.graded", param="numeric"),
          function(m, param) {
            param[2:m@numOutcomes]
          })

setMethod("rpf.setLocation", signature(m="rpf.1dim.graded", param="numeric", loc="numeric"),
          function(m, param, loc) {
              param[2:m@numOutcomes] <- loc
              param
          })

### mdim

setMethod("rpf.rparam", signature(m="rpf.mdim.graded"),
          function(m) {
              a <- rlnorm(m@dimensions, meanlog=m@a.prior.meanlog,
                          sdlog=m@a.prior.sdlog)
              b <- sort(rnorm(m@numOutcomes-1))
              c(a=a,b=b)
          })

setMethod("rpf.startingParam", signature(m="rpf.mdim.graded"),
          function(m) {
            c(a=rep(1,m@dimensions), b=rep(0, (m@numOutcomes-1)))
          })

setMethod("rpf.getLocation", signature(m="rpf.mdim.graded", param="numeric"),
          function(m, param) {
            param[-1:-m@dimensions]
          })

setMethod("rpf.setLocation", signature(m="rpf.mdim.graded", param="numeric", loc="numeric"),
          function(m, param, loc) {
            param[-1:-m@dimensions] <- loc
            param
          })

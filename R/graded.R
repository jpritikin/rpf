### 1dim

setMethod("rpf.rparam", signature(m="rpf.1dim.graded"),
          function(m) {
              a <- rlnorm(1, meanlog=0, sdlog=.5)
              b <- sort(rnorm(m@outcomes-1))
              c(a=a,b=b)
          })

### mdim

setMethod("rpf.rparam", signature(m="rpf.mdim.graded"),
          function(m) {
              a <- rlnorm(m@factors, meanlog=0, sdlog=.5)
              b <- rnorm(m@outcomes-1)
              b <- b[order(-b)]
              c(a=a,b=b)
          })

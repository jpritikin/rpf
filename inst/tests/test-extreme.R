library(rpf)
library(testthat)

spec <- list()
param <- list()
spec [[length(spec) +1]] <- rpf.drm(poor=TRUE)
param[[length(param)+1]] <- c(1, 0, 0)
spec [[length(spec) +1]] <- rpf.drm(multidimensional=TRUE)
param[[length(param)+1]] <- c(1, 0, 0)
spec [[length(spec) +1]] <- rpf.gpcm(3)
param[[length(param)+1]] <- c(1, -10, 10)

# To debug, set breakpoint on Rf_error 
where <- seq(-1000, 1000, 100)
for (ix in 1:length(spec)) {
  ispec <- spec[[ix]]
  iparam <- param[[ix]]
  for (wh in where) rpf.prob(ispec, iparam, wh)
  for (wh in where) rpf.logprob(ispec, iparam, wh)
  w <- rchisq(ispec@numOutcomes, df=6)
  for (wh in where) rpf.gradient(ispec, iparam, wh, w)
}

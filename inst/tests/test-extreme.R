library(rpf)
library(testthat)

context("extremes")

spec <- list()
param <- list()
# repair the poor version of drm TODO
#spec [[length(spec) +1]] <- rpf.drm(poor=TRUE)
#param[[length(param)+1]] <- c(1, 0, 0)
spec [[length(spec) +1]] <- rpf.drm(multidimensional=TRUE)
param[[length(param)+1]] <- c(1, 0, .05, .95)

spec [[length(spec) +1]] <- rpf.grm(3, multidimensional=TRUE)
param[[length(param)+1]] <- c(1, 1, -1)

spec [[length(spec) +1]] <- rpf.nrm(3)
param[[length(param)+1]] <- c(1,  .5, .6, 0, -.6)

# To debug, set breakpoint on Rf_error 
where <- seq(-1000, 1000, 100)
for (ix in 1:length(spec)) {
  ispec <- spec[[ix]]
  iparam <- param[[ix]]
  test_that(paste("extreme values in", class(ispec)), {
    for (wh in where) {
      v <- rpf.prob(ispec, iparam, wh)
      expect_equal(sum(v), 1)
    }
    for (wh in where) {
      v <- rpf.logprob(ispec, iparam, wh)
      expect_equal(sum(exp(v)), 1)
    }
    w <- rchisq(ispec@outcomes, df=6)
    for (wh in where) rpf.dLL(ispec, iparam, wh, w)
  })
}

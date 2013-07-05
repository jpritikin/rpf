library(testthat)
library(rpf)
library(mirt)

context("outfit/infit")

test_that("mirt", {
  set.seed(1)
  a <- matrix(rlnorm(20, meanlog=0, sdlog = .1),ncol=1)
  d <- matrix(rnorm(20),ncol=1)
  data <- simdata(a,d, 1000, rep('dich', 20))
  raschfit <- mirt(data, 1, itemtype='Rasch', D=1, verbose=FALSE)
  coef(raschfit)  # item parameters
  mirt.fit <- itemfit(raschfit)
  scores.full <- fscores(raschfit, full.scores=TRUE)

  spec <- list()
  spec[1:20] <- rpf.drm(multidimensional=FALSE)
  params <- t(simplify2array(coef(raschfit)[1:20]))[,1:4]
  params[,2] <- params[,2] / -params[,1]
  scores <- scores.full[,'F1']
  data.f <- as.data.frame(lapply(as.data.frame(data), ordered))
  fit <- rpf.1dim.fit(spec, t(params), data.f, scores, 2)

  expect_equal(mirt.fit$infit, fit$infit, tolerance=10^-4)
  expect_equal(mirt.fit$outfit, fit$outfit, tolerance=10^-4)
  expect_equal(mirt.fit$z.infit, fit$infit.z, tolerance=10^-3)
                                        # z.outfit TODO
})

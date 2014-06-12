#options(error = browser)
library(testthat)
library(rpf)

context("sumscore")

test_that("observedSumScore", {
  set.seed(1)
  spec <- list()
  spec[1:3] <- rpf.grm(outcomes=3)
  param <- sapply(spec, rpf.rparam)
  data <- rpf.sample(5, spec, param)
  colnames(param) <- colnames(data)
  grp <- list(spec=spec, param=param, data=data)
  obs <- observedSumScore(grp, rep(TRUE, length(spec)))
  expect_equal(obs, c(1L, 1L, 0L, 1L, 1L, 0L, 1L))
  
  dperm <- sample.int(3)
  data <- data[,dperm]
  
  mask <- c(TRUE, FALSE, TRUE)
  obs <- observedSumScore(grp, mask)
  expect_equal(obs, rep(1L, 5))
})

test_that("itemOutcomeBySumScore", {
  set.seed(1)
  spec <- list()
  spec[1:3] <- rpf.grm(outcomes=3)
  param <- sapply(spec, rpf.rparam)
  data <- rpf.sample(5, spec, param)
  colnames(param) <- colnames(data)
  grp <- list(spec=spec, param=param, data=data)
  tbl <- itemOutcomeBySumScore(grp, c(FALSE,TRUE,TRUE), 1L)
  
  want <- structure(c(1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L,  0L, 1L),
                    .Dim = c(5L, 3L), .Dimnames = list(c("0", "1", "2",  "3", "4"), NULL))
  expect_equal(tbl, want)  
})

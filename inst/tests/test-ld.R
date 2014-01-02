library(testthat)
library(rpf)
library(mvtnorm)

context("chen & thissen 1997")

test_that("ptw", {
  set.seed(1)
  
  spec <- list()
  spec[1:6] <- rpf.grm(factors=2)
  gen.param <- sapply(spec, rpf.rparam)
  colnames(gen.param) <- paste("i", 1:ncol(gen.param), sep="")
  gen.param[2,] <- c(0,0,.5,.5,1,1)
  
  theta <- rmvnorm(1000, c(0,0), diag(2))
  resp <- rpf.sample(t(theta), spec, gen.param)
  
  # hide latent factor that we don't know about
  tspec <- list()
  tspec[1:length(spec)] <- rpf.grm(factors=1)
  
  grp <- list(spec=tspec, param=gen.param[-2,], mean=c(0), cov=diag(1), data=resp)
  
  got <- chen.thissen.1997(grp)
  #cat(deparse(round(got$pval[!is.na(got$pval)], 2)))
  expect_equal(got$pval[!is.na(got$pval)], 
               c(-0.33, -0.24, -1.19, -3.52, -0.54, -0.17, 0.21, -3.01, 1.15,
                 0.43, 3.27, 0.73, 6.13, 4.19, 8.41), .001)
  #cat(deparse(round(got$gamma[!is.na(got$gamma)], 3)))
  expect_equal(got$gamma[!is.na(got$gamma)],
               c(-0.064, -0.056, -0.099, -0.002, -0.078, -0.022, 0.04, -0.002,
                 0.061, 0.021, 0.085, 0.043, 0.18, 0.168, 0.3), .01)
})

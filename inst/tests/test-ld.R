library(testthat)
library(rpf)
library(mvtnorm)

context("chen & thissen 1997")

test_that("ct1997", {
  set.seed(1)
  
  spec <- list()
  spec[1:6] <- rpf.grm(factors=2)
  gen.param <- sapply(spec, rpf.rparam, version=1)
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
               c(-1.52, -1.74, -2.54, -6.01, -2.74, -1.06, 1.47, -5.77, 2.42,  1.77, 7.01,
                 2.24, 11.01, 6.99, 21.95), .001)
  #cat(deparse(round(got$gamma[!is.na(got$gamma)], 3)))
  expect_equal(got$gamma[!is.na(got$gamma)],
               c(-0.064, -0.056, -0.099, -0.002, -0.078, -0.022, 0.04, -0.002,
                 0.061, 0.021, 0.085, 0.043, 0.18, 0.168, 0.3), .01)
})

drawRandomProportion <- function(expected) {
  total <- sum(expected)
  prob <- expected / total
  sim <- rep(NA, length(expected))
  rowSim <- sample.int(length(expected), size=total, prob=prob, replace=TRUE)
  sim <- tabulate(rowSim, length(expected))
  sim
}

if (0) {
  spec <- list()
  spec[1:6] <- rpf.grm(factors=1)
  gen.param <- sapply(spec, rpf.rparam)
  grp <- list(spec=spec, param=gen.param, mean=c(0), cov=diag(1))

  pair <- c(1L,2L)
  numPeople <- 185
  E <- (numPeople * pairwiseExpected(grp, pair))

  ms <- pairwiseItemDistribution(grp, pair)
  hist(ms)
  quantile(ms, c(.95, .99))   # 25 38
  
  trials <- 10000
  ms <- rep(NA, trials)
  for (tx in 1:trials) {
    O <- drawPairwiseSample(grp, pair, numPeople)
    ms[tx] <- sum((c(E) - O)^2)
  }

  hist(ms)
  quantile(ms, c(.95, .99))   # 25 38
}

pearson.gof <- function(observed, expected, df) {
  x2 <- sum((observed - expected)^2/expected)
  if (missing(df)) {
    df <- (dim(observed)[1]-1) * (dim(observed)[2]-1)
  }
  pchisq(x2, df, lower.tail=FALSE)
}

if (0) {
  spec <- list()
  spec[1:6] <- rpf.grm(factors=1)
  gen.param <- sapply(spec, rpf.rparam)
  grp <- list(spec=spec, param=gen.param, mean=c(0), cov=diag(1))

  pair <- c(1L,2L)
  numPeople <- 500
  E <- (numPeople * pairwiseExpected(grp, pair))

  trials <- 50
  got <- expand.grid(trial=1:trials, method=c("pearson","rms"), pval=NA)
  for (rep in 1:trials) {
    O <- drawPairwiseSample(grp, pair, 50, qpts=11L, qwidth=4)  # doesn't work!
    O <- O * numPeople / 50
    got[got$trial==rep,'pval'] <- c(pearson.gof(O, E),
                                    pairwiseItemTest(grp, pair, O, qpts=11, qwidth=4))
#                                    ptw2011.gof.test(O, E)
   print(rep)
  }
  #
  require(ggplot2)
 mask <- got$method=="rms"
  #mask <- rep(TRUE, nrow(got))
 pval <- got[mask,'pval']
 got[mask,'pval'] <- 1 / (1+exp(-(logit(pval) - 2.8)))
  
  tbl <- expand.grid(alpha=seq(0,1,.01), method=c("pearson","rms"), emp=NA)
  tbl[tbl$method=="rms",'emp'] <-
    Vectorize(function(x) sum(got[got$method=="rms",'pval'] < x))(seq(0,1,.01))
  tbl[tbl$method=="pearson",'emp'] <-
    Vectorize(function(x) sum(got[got$method=="pearson",'pval'] < x))(seq(0,1,.01))
  tbl$emp <- tbl$emp / trials
  
  ggplot(tbl, aes(alpha, emp, color=method)) + geom_line() +
    geom_abline(slope=1, color="yellow")+ coord_fixed()
}


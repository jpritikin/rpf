library(testthat)
library(plink)
library(rpf)

context("plink")

#source("rpf/R/classes.R")
#source("rpf/R/drm.R")
#source("rpf/R/gpcm.R")

set.seed(1)

theta <- rnorm(2)
theta.2d <- array(dim=c(2,3), data=rnorm(6))

checkDim <- function(item, param) {
    expect_equal(length(param), rpf.numParam(item), class(item))
    expect_equal(length(rpf.rparam(item)), rpf.numParam(item), class(item))
}

test_that("3PL ICC", {
  i1 <- rpf.drm()
  i1.p <- rpf.rparam(i1)
  i1.p[-1] <- 1  # plink doesn't do upper bounds
  checkDim(i1,i1.p)
  expect_equal(drm(t(i1.p[1:3]), theta, 1)@prob[,2],
               rpf.prob(i1, i1.p, theta)[,2],
               "3PL")

  m1 <- rpf.drm(factors=2)
  m1.p <- rpf.rparam(m1)
  m1.p[-1] <- 1  # plink doesn't do upper bounds
  checkDim(m1,m1.p)
  expect_equivalent(rpf.prob(m1, m1.p, theta.2d)[2,],
                    drm(t(m1.p[1:4]), t(theta.2d), dimensions=2)@prob[,3],
                    "M3PL")
})

# Rewrite in terms of nominal model TODO
#
## test_that("GPCM ICC", {
##   i2 <- rpf.gpcm(outcomes=3)
##   i2.p <- rpf.rparam(i2)
##   checkDim(i2,i2.p)
##   expect_equivalent(t(gpcm(t(i2.p), i2@outcomes, theta)@prob[,-1]),
##                     rpf.prob(i2, i2.p, theta),
##                     "GPCM")

# Rewrite in terms of nominal model TODO
#
# i3 <- rpf.gpcm(factors=2, outcomes=3)
# i3.p <- rpf.rparam(i3)
# checkDim(i3,i3.p)
# expect_equivalent(rpf.prob(i3, i3.p, theta.2d),
#     as.matrix(gpcm(t(i3.p),factors=2,cat=3,theta.2d)@prob[,-1:-2]),
#                    "M-GPCM")
#})

# broken, different parameterization TODO
if (0) {
i4 <- rpf.nrm(outcomes=3,factors=2)
i4.p <- rpf.rparam(i4)
#i4.plink <- t(c(i4.p[1:2],ak0=0,i4.p[3:4],g0=0,i4.p[5:6]))
checkDim(i4,i4.p)
expect_equivalent(rpf.prob(i4, i4.p, theta.2d),
                   as.matrix(nrm(x=i4.plink,cat=3, factors=2, t(theta.2d))@prob[,-1:-2]),
                   "NRM")
}

# Not implemented
#
# i5 <- rpf.mcm(outcomes=4,factors=2)
# i5.p <- rpf.rparam(i5)
# checkDim(i5,i5.p)
# expect_equivalent(rpf.prob(i5, i5.p, theta.2d),
#     as.matrix(mcm(t(i5.p),factors=2,cat=4,theta=theta.2d)@prob[,-1:-2]),
#                    "MCM")

test_that("GRM ICC", {
  i6 <- rpf.grm(outcomes=4, multidimensional=FALSE)
  i6.p <- rpf.rparam(i6)
  checkDim(i6,i6.p)
  expect_equivalent(rpf.prob(i6, i6.p, theta),
                    t(grm(t(i6.p), factors=1,cat=4, theta=theta, catprob=TRUE)@prob[,-1]))

  i7 <- rpf.grm(factors=2, outcomes=3)
  i7.p <- rpf.rparam(i7)
  checkDim(i7,i7.p)
  expect_equivalent(rpf.prob(i7, i7.p, theta.2d),
                    t(grm(t(i7.p),dimensions=2,cat=3,t(theta.2d),catprob=TRUE)@prob[,-1:-2]),
                    "M-GRM")
})

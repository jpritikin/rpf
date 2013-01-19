library(testthat)
library(plink)
library(rpf)

#source("rpf/R/classes.R")
#source("rpf/R/drm.R")
#source("rpf/R/gpcm.R")

set.seed(1)

theta <- rnorm(2)
theta.2d <- array(dim=c(3,2), data=rnorm(6))

checkDim <- function(item, param) {
    expect_equal(length(param), item@numParam, class(item))
    expect_equal(length(rpf.rparam(item)), item@numParam, class(item))
    expect_equal(item@numParam, rpf.numParam(item))
}

i1 <- rpf.drm()
i1.p <- rpf.rparam(i1)
checkDim(i1,i1.p)
expect_equal(drm(t(i1.p), theta, 1)@prob[,2],
                   rpf.prob(i1, i1.p, theta)[,2],
                   "3PL")

m1 <- rpf.drm(dimensions=2)
m1.p <- rpf.rparam(m1)
checkDim(m1,m1.p)
expect_equivalent(rpf.prob(m1, m1.p, theta.2d)[,2],
    drm(t(m1.p), theta.2d, dimensions=2)@prob[,3],
                   "M3PL")

i2 <- rpf.gpcm(numOutcomes=3)
i2.p <- rpf.rparam(i2)
checkDim(i2,i2.p)
expect_equivalent(as.matrix(gpcm(t(i2.p), i2@numOutcomes, theta)@prob[,-1]),
                   rpf.prob(i2, i2.p, theta),
                   "GPCM")

i3 <- rpf.gpcm(dimensions=2, numOutcomes=3)
i3.p <- rpf.rparam(i3)
checkDim(i3,i3.p)
expect_equivalent(rpf.prob(i3, i3.p, theta.2d),
    as.matrix(gpcm(t(i3.p),dimensions=2,cat=3,theta.2d)@prob[,-1:-2]),
                   "M-GPCM")

i4 <- rpf.nrm(numOutcomes=3,dimensions=2)
i4.p <- rpf.rparam(i4)
checkDim(i4,i4.p)
expect_equivalent(rpf.prob(i4, i4.p, theta.2d),
                   as.matrix(nrm(x=t(i4.p),cat=3, dimensions=2, theta.2d)@prob[,-1:-2]),
                   "NRM")

i5 <- rpf.mcm(numOutcomes=4,dimensions=2)
i5.p <- rpf.rparam(i5)
checkDim(i5,i5.p)
expect_equivalent(rpf.prob(i5, i5.p, theta.2d),
    as.matrix(mcm(t(i5.p),dimensions=2,cat=4,theta=theta.2d)@prob[,-1:-2]),
                   "MCM")

i6 <- rpf.grm(numOutcomes=4)
i6.p <- rpf.rparam(i6)
checkDim(i6,i6.p)
expect_equivalent(rpf.prob(i6, i6.p, theta),
  as.matrix(grm(t(i6.p), dimensions=1,cat=4, theta=theta, catprob=TRUE)@prob[,-1]))

i7 <- rpf.grm(dimensions=2, numOutcomes=3)
i7.p <- rpf.rparam(i7)
checkDim(i7,i7.p)
expect_equivalent(rpf.prob(i7, i7.p, theta.2d),
                   as.matrix(grm(t(i7.p),dimensions=2,cat=3,theta.2d,catprob=TRUE)@prob[,-1:-2]),
                   "M-GRM")

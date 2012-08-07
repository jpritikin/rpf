library(RUnit)
library(plink)
library(rpf)

#source("rpf/R/classes.R")
#source("rpf/R/drm.R")
#source("rpf/R/gpcm.R")

set.seed(1)

theta <- rnorm(2)
theta.2d <- array(dim=c(2,2), data=rnorm(4))

i1 <- rpf.drm()
i1.p <- rpf.rparam(i1)
checkEqualsNumeric(drm(i1.p, theta, 1)@prob[,2],
                   rpf.prob(i1, i1.p, theta)[,2],
                   "3PL")

m1 <- rpf.drm(dimensions=2)
m1.p <- rpf.rparam(m1)
checkEqualsNumeric(rpf.prob(m1, m1.p, theta.2d)[,2],
    drm(m1.p, theta.2d, dimensions=2)@prob[,3],
                   "M3PL")

i2 <- rpf.gpcm(numOutcomes=3)
i2.p <- rpf.rparam(i2)
checkEqualsNumeric(as.matrix(gpcm(i2.p, i2@numOutcomes, theta)@prob[,-1]),
                   rpf.prob(i2, i2.p, theta),
                   "GPCM")

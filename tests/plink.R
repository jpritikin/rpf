library(RUnit)
library(plink)
library(rpf)

#source("rpf/R/classes.R")
#source("rpf/R/drm.R")
#source("rpf/R/gpcm.R")

set.seed(1)

theta <- rnorm(2)

i1 <- rpf.drm()
i1.p <- rpf.rparam(i1)
checkEqualsNumeric(drm(i1.p, theta, 1)@prob[,2],
                   rpf.prob(i1, c(i1.p), theta)[,2],
                   "3PL")

i2 <- rpf.gpcm(numOutcomes=3)
i2.p <- rpf.rparam(i2)
checkEqualsNumeric(as.matrix(gpcm(i2.p, i2@numOutcomes, theta)@prob[,-1]),
                   rpf.prob(i2, c(i2.p), theta),
                   "GPCM")

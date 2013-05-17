library(testthat)
library(rpf)
library(mirt)
#options(error = utils::recover)

i.count <- 5
spec <- list()

spec[1:i.count] <- rpf.drm(multidimensional=TRUE)
data <- rpf.sample(100, spec)
data <- simplify2array(lapply(data, unclass)) - 1

suppressWarnings(fit <- mirt(data, 1, rep('3PL',i.count), D=1, technical=list(NCYCLES=1)))

for (ix in 1:i.count) {
  ii <- extract.item(fit, 1)
  expect_equal(c(probtrace(ii, c(-1,0,1))),
               c(rpf.prob(spec[[1]], ii@par[1:3], c(-1,0,1))))
}

spec[1:i.count] <- rpf.grm(numOutcomes=3, multidimensional=TRUE)

data <- rpf.sample(100, spec)
data <- simplify2array(lapply(data, unclass)) - 1

suppressWarnings(fit <- mirt(data, 1, rep('graded',i.count), D=1, technical=list(NCYCLES=1)))

for (ix in 1:i.count) {
  ii <- extract.item(fit, 1)
  expect_equal(c(probtrace(ii, c(-1,0,1))),
               c(rpf.prob(spec[[1]], ii@par[1:3], c(-1,0,1))))
}

#suppressWarnings(fit <- mirt(data, 1, rep('nominal',i.count), D=1, technical=list(NCYCLES=1)))
#ii <- extract.item(fit, 1)
#str(ii)
#ii@par

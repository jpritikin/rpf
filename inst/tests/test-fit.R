#options(error = browser)
library(testthat)
library(rpf)

context("ot2000")

# Thissen, Pommerich, Billeaud, & Williams (1995)
test_that("tpbw1995-table2", {
  spec <- list()
  spec[1:3] <- rpf.grm(outcomes=4)
  
  param <- matrix(c(1.87, .65, 1.97, 3.14,
                    2.66, .12, 1.57, 2.69,
                    1.24, .08, 2.03, 4.3), nrow=4)
  # fix parameterization
  param <- apply(param, 2, function(p) c(p[1], p[2:4] * -p[1]))
  
  grp <- list(spec=spec, mean=0, cov=matrix(1,1,1), param=param)
  
  got <- sumScoreEAP(grp)
  
  expect_equal(sum(got[,1]), 1, tolerance=.001)
  
  #cat(deparse(round(got[,2],3)))
  ssP <- c(0.325, 0.241, 0.183, 0.123, 0.069, 0.035, 0.016, 0.006, 0.002,  0)
  expect_equal(got[,1], ssP, tolerance=.01)
  ssEAP <- c(-0.885, -0.179, 0.332, 0.744, 1.115, 1.482, 1.843, 2.212, 2.622,  2.999)
  expect_equal(got[,2], ssEAP, tolerance=.01)
  ssVar <- c(0.494, 0.378, 0.329, 0.299, 0.297, 0.296, 0.29, 0.296, 0.313,  0.328)
  expect_equal(got[,3], ssVar, tolerance=.01)
})

verifySumP <- function(grp, sseap, N=2000) {  # a good fit is close to 1
  sim <- apply(sapply(rpf.sample(N, grp=grp), unclass), 1, function(r) sum(r-1))
  observed <- tabulate(1+sim, length(sseap[,1]))
#  print(observed/N)
  ptw2011.gof.test(observed, N*sseap[,1])
}

if (0) {
  fm <- read.flexmirt("~/ifa/ifa-2d-mg/2d-mg-prm.txt")
  
  got <- sumScoreEAP(fm$G1, 5, 21L)  # matches flexmirt exactly
  verifySumP(fm$G1, got)
  
  got <- sumScoreEAP(fm$G2, 5, 21L)  # doesn't match flexmirt
  verifySumP(fm$G2, got, N=5000)  # but looks feasible
  
  got <- sumScoreEAP(fm$G3, 5, 21L)  # doesn't match flexmirt
  verifySumP(fm$G3, got, N=5000)  # but looks feasible
}

if (0) {
  # cai2009
  fm <- structure(list(G1 = structure(list(param = structure(c(0.992675,  0.646717, 0, 0, 0.876469, 1.41764, 1.25402, 0, 0, 0.0826927,  1.76547, 1.20309, 0, 0, -0.346706, 2.1951, 0.844399, 0, 0, -0.978301,  1.37774, 0, 1.06694, 0, 0.992373, 1.80365, 0, 0.814109, 0, 0.213559,  2.15718, 0, 1.58086, 0, -0.418129, 1.18201, 0, 1.56533, 0, -1.24173,  1.80474, 0, 0, 1.0774, 0.810718, 2.60754, 0, 0, 1.23507, 0.0598008,  1.01874, 0, 0, 0.724402, -0.294029, 1.68916, 0, 0, 1.37546, -1.13333 ), .Dim = c(5L, 12L), .Dimnames = list(NULL, c("i1", "i2", "i3",  "i4", "i5", "i6", "i7", "i8", "i9", "i10", "i11", "i12"))), mean = structure(c(0.822622,  -0.290462, 0.19672, 0.733993), .Names = c("X6", "X7", "X8", "X9" )), cov = structure(c(0.826046, 0, 0, 0, 0, 1.656, 0, 0, 0, 0,  1.11263, 0, 0, 0, 0, 1.07878), .Dim = c(4L, 4L))), .Names = c("param",  "mean", "cov")),
                       G2 = structure(list(param = structure(c(0.992675,  0.646717, 0, 0, 0, 0.876469, 1.41764, 1.25402, 0, 0, 0, 0.0826927,  1.76547, 1.20309, 0, 0, 0, -0.346706, 2.1951, 0.844399, 0, 0,  0, -0.978301, 1.37774, 0, 1.06694, 0, 0, 0.992373, 1.80365, 0,  0.814109, 0, 0, 0.213559, 2.15718, 0, 1.58086, 0, 0, -0.418129,  1.18201, 0, 1.56533, 0, 0, -1.24173, 1.80474, 0, 0, 1.0774, 0,  0.810718, 2.60754, 0, 0, 1.23507, 0, 0.0598008, 1.01874, 0, 0,  0.724402, 0, -0.294029, 1.68916, 0, 0, 1.37546, 0, -1.13333,  1.75531, 0, 0, 0, 1.20652, 0.875564, 1.26308, 0, 0, 0, 1.25013,  0.196607, 1.44526, 0, 0, 0, 0.990354, -0.351181, 1.89461, 0,  0, 0, 0.85611, -1.09382), .Dim = c(6L, 16L), .Dimnames = list(     NULL, c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", "i9",      "i10", "i11", "i12", "i13", "i14", "i15", "i16"))), mean = structure(c(0,  0, 0, 0, 0), .Names = c("X6", "X7", "X8", "X9", "X10")), cov = structure(c(1,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 1), .Dim = c(5L, 5L))), .Names = c("param", "mean", "cov" ))), .Names = c("G1", "G2"))
  spec <- list()
  spec[1:ncol(fm$G2$param)] <- rpf.grm(factors = 5)
  fm$G2$spec <- spec
  got <- sumScoreEAP(fm$G2, 5, 21L)
  verifySumP(fm$G2, got, N=3000)  # wrong
}

test_that("simple case", {
  set.seed(1)

  spec <- list()
  spec[1:3] <- rpf.drm()

  gen.p <- matrix(c(1,0,0,1,
                    1,-1,0,1,
                    1,1,0,1), ncol=3, nrow=4)  # note byrow=FALSE
  gen.p[3:4,] <- logit(gen.p[3:4,])
  data <- rpf.sample(200, spec, gen.p)
                                        #  write.csv(data, "fit-test.csv", row.names=FALSE, quote=FALSE)

  param <- matrix(c(.71,.03,0, 1,
                    1.53,-.85, 0, 1,
                    1.1,.9,0, 1), ncol=3, nrow=4)
  param[3:4,] <- logit(param[3:4,])
  
  colnames(param) <- paste("i", 1:3, sep="")
  grp <- list(spec=spec, param=param, data=data, mean=0, cov=matrix(1,1,1))

  got <- rpf.SitemFit(grp, method="pearson")
  
  tbl <- round(t(sapply(got, function(row) c(stat=row$statistic, df=row$df, p=row$pval))),2)
  expect_equal(sum(tbl[,'df']), 9)  # not sure TODO
  expect_equal(sum(tbl[,'stat']), 22.55)
})

test_that("orlando-thissen-2000", {
  require(rpf)
  set.seed(7)
  grp <- list(spec=list())
  grp$spec[1:20] <- rpf.grm()
  grp$param <- sapply(grp$spec, rpf.rparam)
  colnames(grp$param) <- paste("i", 1:20, sep="")
  grp$mean <- 0
  grp$cov <- diag(1)
  grp$free <- grp$param != 0
  grp$data <- rpf.sample(500, grp=grp)
  
  got <- rpf.SitemFit(grp, method="pearson")
  
  E1orig <- structure(c(2, 3.99, 12.95, 21.82, 21.66, 25.29, 30.5, 32.29,  35.09, 24.41, 34.33, 15.74,
              18.61, 14.44, 10.79, 5.72, 3.82,  1.5, 0.23, 0.12, 0, 0.01, 0.05, 0.18, 0.34,
              0.71, 1.5, 2.71,  4.91, 5.59, 12.67, 9.26, 17.39, 21.56, 26.21, 23.28, 27.18,
              19.5,  5.77, 5.88), .Dim = c(20L, 2L))
  expect_equal(got[[1]]$orig.expected, E1orig, tolerance=.01)
  
  E1 <- structure(c(2, 3.99, 12.95, 21.82, 21.66, 25.29, 30.5, 32.29,  35.09, 24.41, 34.33,
                    15.74, 18.61, 14.44, 10.79, 5.72, 3.82,  1.85, 0, 0, 0, 0, 0, 0, 0, 1.3,
                    1.5, 2.71, 4.91, 5.59, 12.67,  9.26, 17.39, 21.56, 26.21, 23.28, 27.18,
                    19.5, 5.77, 5.88), .Dim = c(20L,  2L))
  mask <- !is.na(got[[1]]$expected)
  expect_equal(got[[1]]$expected[mask], E1[mask], tolerance=.01)
  
  expect_equal(got[[1]]$pval, 0.3931, tolerance=.001)  # not sure about df TODO
})

if (0) {
  library(mirt)
  dat <- expand.table(LSAT6)
  model <- mirt.model('F = 1-5
                    CONSTRAIN = (1-5, a1)')
  (mod <- mirt(dat, model))
  itemfit(mod, X2=TRUE, S_X2.tables=TRUE)$E[[3]]
  itemfit(mod, X2=TRUE)
  
  library(rpf)
  spec <- list()
  spec[1:5] <- rpf.drm()
  param <- sapply(coef(mod)[-6], function(x) x)
  param[3:4,] <- rpf::logit(param[3:4,])
  dat2 <- as.data.frame(lapply(as.data.frame(dat), ordered, levels=0:1))
  grp <- list(spec=spec, param=param, mean=0, cov=diag(1), data=dat2)
  got <- rpf.SitemFit(grp, method="pearson", alt=FALSE)
}

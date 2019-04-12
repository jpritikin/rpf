suppressWarnings(RNGversion("3.5"))
#options(error = browser)
library(testthat)
library(rpf)

context("ot2000")

test_that("simple case", {
  require(rpf)
  set.seed(1)

  spec <- list()
  spec[1:3] <- list(rpf.drm())

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

  got <- SitemFit(grp, method="pearson")
  
  tbl <- round(t(sapply(got, function(row) c(stat=row$statistic, df=row$df, p=row$pval))),2)
  expect_equal(sum(tbl[,'df']), 9)  # not sure TODO
  expect_equal(sum(tbl[,'stat']), 22.55)
})

test_that("orlando-thissen-2000", {
  require(rpf)
  set.seed(7)
  grp <- list(spec=list())
  grp$spec[1:20] <- list(rpf.grm())
  grp$param <- sapply(grp$spec, rpf.rparam, version=1L)
  colnames(grp$param) <- paste("i", 1:20, sep="")
  grp$mean <- 0
  grp$cov <- diag(1)
  grp$free <- grp$param != 0
  grp$data <- rpf.sample(500, grp=grp)
  
  got <- SitemFit(grp, method="pearson", omit=1)
  expect_true(all(sapply(got, function(ii) is.null(ii$omitted))))
  stat <- sapply(got, function(x) x$statistic)
  names(stat) <- NULL
  Estat <- c(28.44, 14.32, 20.12, 16.28, 15.16, 10.27, 41.52, 9.33, 9.22, 
             16.99, 13.94, 11.99, 7.19, 10.47, 11.76, 13.53, 22.21, 19.72,
             14.82, 19.92)
  expect_equal(stat, Estat, tolerance=.01)
  
  E1orig <- structure(c(2, 3.99, 12.95, 21.82, 21.66, 25.29, 30.5, 32.29,  35.09, 24.41, 34.33, 15.74,
              18.61, 14.44, 10.79, 5.72, 3.82,  1.5, 0.23, 0.12, 0, 0.01, 0.05, 0.18, 0.34,
              0.71, 1.5, 2.71,  4.91, 5.59, 12.67, 9.26, 17.39, 21.56, 26.21, 23.28, 27.18,
              19.5,  5.77, 5.88), .Dim = c(20L, 2L))
  oexp <- got[[1]]$orig.expected
  dimnames(oexp) <- NULL
  expect_equal(oexp, E1orig, tolerance=.01)
  
  E1 <- structure(c(2, 3.99, 12.95, 21.82, 21.66, 25.29, 30.5, 32.29,  35.09, 24.41, 34.33,
                    15.74, 18.61, 14.44, 10.79, 5.72, 3.82,  1.85, 0, 0, 0, 0, 0, 0, 0, 1.3,
                    1.5, 2.71, 4.91, 5.59, 12.67,  9.26, 17.39, 21.56, 26.21, 23.28, 27.18,
                    19.5, 5.77, 5.88), .Dim = c(20L,  2L))
  mask <- !is.na(got[[1]]$expected)
  expect_equal(got[[1]]$expected[mask], E1[mask], tolerance=.01)
  
  expect_equal(got[[1]]$pval, -2.88, tolerance=.01)

  got <- SitemFit(grp, method="pearson", alt=TRUE)
  stat <- sapply(got, function(x) x$statistic)
  names(stat) <- NULL
  #cat(deparse(round(stat, 2)))
  Estat <- c(16.6, 13.68, 13.19, 14.03, 19.86, 11.02, 29.92, 6.95, 18.31,
             14.78, 10.06, 9.27, 5.67, 9.3, 14.75, 18.03, 20.18, 20.19, 15.89,  11.57)
  expect_equal(stat, Estat, tolerance=.01)
})

test_that("fit w/ mcar", {
  require(rpf)
  require(testthat)
  set.seed(7)
  grp <- list(spec=list(), qwidth=5, qpoints=31)
  grp$spec[1:20] <- list(rpf.grm())
  grp$param <- sapply(grp$spec, rpf.rparam, version=1L)
  colnames(grp$param) <- paste("i", 1:20, sep="")
  grp$free <- grp$param != 0
  grp$data <- rpf.sample(500, grp=grp, mcar=.1)

  got <- sumScoreEAPTest(omitMostMissing(grp, 3L))
  expect_equal(got$n, 101L)
  expect_equal(got$rms.p, -1.87, tolerance=.01)
  expect_equal(got$pearson.df, 17L)
  expect_equal(got$pearson.p, -1.26, tolerance=.01)

  grp1 <- grp
  grp1$data <- grp$data[1:250,]
  grp2 <- grp
  grp2$data <- grp$data[251:500,]
  e1 <- sumScoreEAPTest(grp1) + sumScoreEAPTest(grp2)
  e2 <- sumScoreEAPTest(grp)
  chk <- c('n','pearson.df', 'pearson.chisq', 'pearson.p')
  expect_equal(unlist(e1[chk]), unlist(e2[chk]))
  
  got <- SitemFit(grp, omit=2L)
  stat <- sapply(got, function(x) x$statistic)
  names(stat) <- NULL
  Estat <- c(14.78, 8.5, 12.88, 15.42, 16.39, 20.44, 17.88, 46.23, 8.32,
             13, 12.61, 9.47, 8.77, 10.35, 12.57, 17.4, 15.18, 14.57,
             16.83,  14.14)
  expect_equal(stat, Estat, tolerance=.01)
  
  e1 <- SitemFit(grp1) + SitemFit(grp2)
  e2 <- SitemFit(grp)
  expect_equal(sapply(e1, function(ii) ii$statistic),
               sapply(e2, function(ii) ii$statistic))
})

test_that("2tier fit", {
  set.seed(1)
  require(rpf)
  numItems <- 6
  spec <- list()
  spec[1:numItems] <- list(rpf.drm(factors=3))
  param <- sapply(spec, rpf.rparam, version=1)
  gsize <- numItems/3
  for (gx in 0:2) {
    if (gx != 1) {
      param['a2', seq(gx * gsize+1, (gx+1)*gsize)] <- 0
    }
    if (gx != 2) {
      param['a3', seq(gx * gsize+1, (gx+1)*gsize)] <- 0
    }
  }
  grp <- list(spec=spec, param=param, mean=runif(3, -1, 1), cov=diag(runif(3,.5,2)))
  grp$data <- rpf.sample(500, grp=grp)
  colnames(grp$param) <- colnames(grp$data)

  got <- SitemFit(grp, qwidth=2, qpoints=5L, .twotier=FALSE)
  tt <- SitemFit(grp, qwidth=2, qpoints=5L, .twotier=TRUE)
  expect_equal(sapply(got, function(x) x$statistic),
               sapply(tt, function(x) x$statistic), .001)

  got <- SitemFit(grp, qwidth=2, qpoints=5L, .twotier=FALSE, alt = TRUE)
  tt <- SitemFit(grp, qwidth=2, qpoints=5L, .twotier=TRUE, alt=TRUE)
  expect_equal(sapply(got, function(x) x$statistic),
               sapply(tt, function(x) x$statistic), .001)
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
  spec[1:5] <- list(rpf.drm())
  param <- sapply(coef(mod)[-6], function(x) x)
  param[3:4,] <- rpf::logit(param[3:4,])
  dat2 <- as.data.frame(lapply(as.data.frame(dat), ordered, levels=0:1))
  grp <- list(spec=spec, param=param, mean=0, cov=diag(1), data=dat2)
  got <- SitemFit(grp, method="pearson", alt=FALSE)
}

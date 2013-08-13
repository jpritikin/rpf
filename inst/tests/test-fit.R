#options(error = utils::recover)
library(testthat)
library(rpf)

context("ot2000")

test_that("simple case", {
  set.seed(1)

  spec <- list()
  spec[1:3] <- rpf.drm()

  gen.p <- matrix(c(1,0,0,1,
                    1,-1,0,1,
                    1,1,0,1), ncol=3, nrow=4)  # note byrow=FALSE
  data <- rpf.sample(200, spec, gen.p)
                                        #  write.csv(data, "fit-test.csv", row.names=FALSE, quote=FALSE)

  param <- matrix(c(.71,.03,0, 1,
                    1.53,-.85, 0, 1,
                    1.1,.9,0, 1), ncol=3, nrow=4)

  got <- rpf.ot2000.chisq(spec, param, param!=0 & param!=1, data)

  tbl <- round(t(sapply(got, function(row) c(chisq=row$statistic, df=row$df, p=row$p.value))),2)
  expect_equal(sum(tbl[,2]), 3)
  expect_equal(sum(tbl[,1]), 22.55)
})

options(error = utils::recover)
library(testthat)
library(rpf)
data(kct)

context("outfit/infit")

responses <- kct.people[,paste("V",2:19, sep="")]
rownames(responses) <- kct.people$NAME
responses <- responses[1:34,4:17]
data <- sapply(responses, unclass) - 1

test_that("kct", {
  scores <- kct.people$MEASURE
  params <- as.data.frame(cbind(1, kct.items$MEASURE, 0))
  rownames(params) <- kct.items$NAME
  items<-list()
  items[1:18] <- rpf.drm(multidimensional=FALSE)

  fit <- rpf.1dim.fit(items[1:14], t(cbind(params[4:17,],1)),
                      responses, scores[1:34], 2)

  expect_equal(fit$infit, kct.items$IN.MSQ[4:17], tolerance=.002)
  expect_equal(fit$infit.z, kct.items$IN.ZSTD[4:17], tolerance=.01)
  expect_equal(fit$outfit, kct.items$OUT.MSQ[4:17], tolerance=.002)

                                        #fit$outfit.z
                                        #fit$outfit.z - kct.items$OUT.ZSTD[4:17]
                                        #expect_equal(fit$outfit.z, kct.items$OUT.ZSTD[4:17])   # TODO
})

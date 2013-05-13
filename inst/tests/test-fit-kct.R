options(error = utils::recover)
library(testthat)
library(rpf)

data(kct)

scores <- kct.people$MEASURE
params <- as.data.frame(cbind(1, kct.items$MEASURE, 0))
rownames(params) <- kct.items$NAME
items<-list()
items[1:18] <- rpf.drm()
responses <- kct.people[,paste0("V",2:19)]
rownames(responses) <- kct.people$NAME

fit <- rpf.1dim.fit(items[1:14], params[4:17,], responses[1:34,4:17], scores[1:34], 2)

expect_equal(fit$infit, kct.items$IN.MSQ[4:17], tolerance=.002)
expect_equal(fit$infit.z, kct.items$IN.ZSTD[4:17], tolerance=.01)
expect_equal(fit$outfit, kct.items$OUT.MSQ[4:17], tolerance=.002)

#fit$outfit.z
#fit$outfit.z - kct.items$OUT.ZSTD[4:17]
#expect_equal(fit$outfit.z, kct.items$OUT.ZSTD[4:17])   # TODO

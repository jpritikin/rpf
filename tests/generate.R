library(RUnit)
library(rpf)

set.seed(1)

### 1 dimensional items

i1 <- rpf.drm()
i1.p <- rpf.rparam(i1)
i2 <- rpf.gpcm(numOutcomes=3)
i2.p <- rpf.rparam(i2)
data <- rpf.sample(3, list(i1,i2), list(i1.p, i2.p))
checkEqualsNumeric(c(data), c(2, 1, 1, 3, 3, 3))

data <- rpf.sample(runif(3), list(i1,i2), list(i1.p, i2.p))
checkEqualsNumeric(c(data), c(2, 1, 2, 2, 2, 2))

### multidimensional items

numItems <- 4
items <- vector("list", numItems)
correct <- vector("list", numItems)

i1 <- rpf.drm(dimensions=2)
i2 <- rpf.drm(dimensions=1, multidimensional=TRUE)

for (ix in 1:(numItems-1)) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1)
}
items[[4]] <- i2
correct[[4]] <- rpf.rparam(i2)

design <- matrix(c(1, 1, 1, 1,
                   2, 2, 3, NA), nrow=2, byrow=TRUE)
data <- rpf.sample(3, items, correct, design)
checkEqualsNumeric(c(data), c(1, 1, 1, 2, 2, 2, 2, 1, 2, 1, 1, 2))

### multidimension, no design

numItems <- 3
items <- vector("list", numItems)
correct <- vector("list", numItems)

i1 <- rpf.drm(dimensions=2)
for (ix in 1:numItems) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1)
}

data <- rpf.sample(3, items, correct)
checkEqualsNumeric(c(data), c(1, 2, 1, 1, 2, 2, 1, 1, 1))

### uneven 2d design

numItems <- 5
numPersons <- 3
maxDim <- 2

items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
	items[[ix]] <- rpf.drm(dimensions=min(ix,maxDim),
                               multidimensional=TRUE)
	correct[[ix]] <- rpf.rparam(items[[ix]])
	correct[[ix]][[4]] <- 0   # no guessing, for now
}

maxParam <- max(vapply(items, function(i) i@numParam, 0))
maxOutcomes <- max(vapply(items, function(i) i@numOutcomes, 0))

design <- matrix(c(1, 1,1,1,2,
		   NA,2,2,2,1), byrow=TRUE, nrow=2)

data <- rpf.sample(numPersons, items, correct, design)

checkEqualsNumeric(c(data),
                   c(1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L))

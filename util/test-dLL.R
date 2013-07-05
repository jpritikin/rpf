library(testthat)
library(OpenMx)  # an unrelease version of OpenMx is required for this test
library(rpf)
library(numDeriv)
options(error = utils::recover)

context("dLL")

unpackHession <- function(deriv, np) {
  hess <- matrix(NA, nrow=np, ncol=np)
  dx <- np+1
  for (hr in 1:np) {
    hess[1:hr,hr] <- hess[hr,1:hr] <- deriv$D[dx:(dx+hr-1)]
    dx <- dx + hr
  }
  hess
}

myseed <- as.integer(runif(1) * 1e7)
#print(paste("set.seed =",myseed))
set.seed(myseed)
#set.seed(3090526)

numItems <- 3
items <- list()
items[[1]] <- rpf.drm(factors=2)
items[[2]] <- rpf.grm(outcomes=3, factors=2)
T.a <- matrix(rnorm(9),3,3) + diag(3)
T.c <- matrix(rnorm(9),3,3) + diag(3)

items[[3]] <- rpf.nrm(outcomes=4, factors=2,
                      T.a=T.a, T.c=T.c)

# If not all outcomes are represented then lots of warnings result.
# This is not a cause for concern.
data <- rpf.sample(200, items)

max.spec.len <- max(vapply(items, function(m) length(m@spec), 0))

spec <- mxMatrix(name="ItemSpec", nrow=max.spec.len, ncol=numItems,
                 values=NA, free=FALSE, byrow=TRUE)
for (ix in 1:length(items)) {
  svec <- items[[ix]]@spec
  spec@values[1:length(svec),ix] <- svec
}

starting <- list(c(1.4, 1, 0, .1, .9),
                 c(1.4, 1, -.5, -1),
                 c(1.4,  1,  rep(0,6)))
starting.len <- max(vapply(starting, length, 0))

ip.mat <- mxMatrix(name="itemParam", nrow=starting.len, ncol=numItems,
                   values=0, free=FALSE)

for (sx in 1:length(starting)) {
  v <- starting[[sx]]
  ip.mat@values[1:length(v),sx] <- v
  ip.mat@free[1:length(v),sx] <- TRUE
}

Eip <- mxMatrix(name="EitemParam", nrow=dim(ip.mat@values)[1], ncol=numItems,
                values=ip.mat@values)

m.mat <- mxMatrix(name="mean", nrow=1, ncol=2, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=2, ncol=2, values=diag(2), free=FALSE)
m2 <- mxModel(model="drm1", ip.mat, spec, Eip, m.mat, cov.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov",
                ItemSpec="ItemSpec",
                EItemParam="EitemParam",
                qpoints=12),
              mxFitFunctionBA81(ItemParam="itemParam", rescale=FALSE),
              mxComputeSequence(steps=list(
                                  mxComputeOnce('expectation', context='E'),
                                  mxComputeOnce('fitfunction', gradient=TRUE, hessian=TRUE)
                                  )))

spoint <- list(c(1.4, 1, 0, .1, .9),
               c(1.4, 1, .5, -.5),
               c(0.89,  0.33,  0.05,  0.18,  0.07, -2.03,  0.17,  0.26))
spoint.len <- vapply(spoint, length, 0)

for (ii in 1:numItems) {
  np <- length(spoint[[ii]])
  m2@matrices$itemParam@values <- Eip@values
  
  m2@matrices$itemParam@values[1:np,ii] <- spoint[[ii]]
  m2 <- mxRun(m2, useOptimizer=FALSE, silent=TRUE)
  
  offset <- 1
  if (ii > 1) offset <- sum(c(1,spoint.len[1:(ii-1)]))
  offset.i <- seq(offset,offset+np-1)

  grad1 <- m2@output$gradient[offset.i]
  names(grad1) <- NULL
  hess <- m2@output$hessian[offset.i, offset.i]
  for (cx in 2:dim(hess)[1]) for (rx in 1:cx) hess[rx,cx] <- hess[cx,rx]
  
  deriv <- genD(function(param) {
    np <- length(param)
    m2@matrices$itemParam@values[1:np,ii] <- param
    m2 <- mxRun(m2, silent=TRUE)
    m2@output$minimum
  }, spoint[[ii]], method.args=list(eps=0.01, d=0.01, r=2))

  emp.hess <- unpackHession(deriv, np)
  emp.hess[is.na(emp.hess)] <- 0

  if (0) {
    print(paste("Item", ii))
    print(grad1)
    print(deriv$D[1:np])
    cat("T.a=",deparse(T.a),"\n")
    cat("T.c=",deparse(T.c),"\n")
    cat("an=",deparse(hess),"\n")
    cat("emp=",deparse(emp.hess),"\n")
    print(round(hess - emp.hess, 2))
  }
  
  expect_equal(deriv$D[1:np], grad1, tolerance=1e-6)
  expect_equal(emp.hess, hess, tolerance=1e-4)
}
#warnings()

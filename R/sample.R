##' Randomly sample response patterns given a list of items
##'
##' Returns a random sample of response patterns given
##' a list of item models and parameters.
##'
##' @name rpf.sample
##' @usage
##' rpf.sample(theta, items, params, design)
##' @param theta either a vector of trait abilities or
##' the number of abilities to draw from N(0,1)
##' @param items a list of item models
##' @param params a list of item parameters. If omitted, random item
##' parameters are generated for each item model.
##' @param design assigns person abilities to item dimensions
##' @return Returns a matrix of response patterns
##' @export
##' @examples
##' # 1 dimensional items
##' i1 <- rpf.drm()
##' i1.p <- rpf.rparam(i1)
##' i2 <- rpf.gpcm(numOutcomes=3)
##' i2.p <- rpf.rparam(i2)
##' rpf.sample(5, list(i1,i2), list(i1.p, i2.p))
##'
##' # multidimensional items
##' numItems <- 4
##' items <- vector("list", numItems)
##' correct <- vector("list", numItems)
##'
##' i1 <- rpf.drm(dimensions=2)
##' i2 <- rpf.drm(dimensions=1, multidimensional=TRUE)
##'
##' for (ix in 1:(numItems-1)) {
##'   items[[ix]] <- i1
##'   correct[[ix]] <- rpf.rparam(i1)
##' }
##' items[[4]] <- i2
##' correct[[4]] <- rpf.rparam(i2)
##' 
##' design <- matrix(c(1, 1, 1, 1,
##'                    2, 2, 3, NA), nrow=2, byrow=TRUE)
##' rpf.sample(10, items, correct, design)
##' @seealso \code{\link{sample}}
rpf.sample <- function(theta, items, params, design) {
  maxDim <- max(vapply(items, function(i) i@dimensions, 0))
  if (missing(design)) {
    if (maxDim > 1) {
      stop("The design matrix must be provided for multidimensional item models")
    }
    design <- 1
  }
  maxAbilities <- max(design, na.rm=TRUE)

  numItems <- length(items)
  if (maxDim > 1 && any(dim(design) != c(maxDim,numItems))) {
    stop(paste("The design matrix must have", maxDim, "rows and ",numItems,"columns"))
  }

  numPeople <- NA
  if (is.numeric(theta) && length(theta) == 1) {
    numPeople <- theta
    theta <- array(rnorm(numPeople * maxAbilities),
                   dim=c(numPeople, maxAbilities))
  } else if (maxDim == 1 && is.vector(theta)) {
    numPeople <- length(theta)
    theta <- array(theta, dim=c(numPeople, maxAbilities))
  } else {
    numPeople <- dim(theta)[1]
  }

  if (missing(params)) {
    params <- lapply(items, rpf.rparam)
  }
  outcomes <- vapply(items, function(i) i@numOutcomes, 0)
  
  ret <- array(dim=c(numPeople, numItems))
  for (ix in 1:numItems) {
    i <- items[[ix]]
    param <- params[[ix]]
    P <- NA
    if (maxDim==1) {
      P <- rpf.prob(i, param, theta)
    } else {
      cols <- design[,ix]
      cols <- cols[!is.na(cols)]
      i.theta <- as.matrix(theta[,cols])
      P <- rpf.prob(i, param, i.theta)
    }
    ret[,ix] <- apply(P, c(1), sample, x=1:i@numOutcomes, size=1, replace=F)
  }
  return(ret)
}

##' Randomly sample response patterns given a list of items
##'
##' Returns a random sample of response patterns given a list of item
##' models and parameters.
##'
##' The design matrix can accomodate more person abilities than item
##' dimension. Refer to Cai (2010) for design matrix examples.
##'
##' TODO: Add restrictions to design matrix to match restrictions
##' imposed by Cai (2010).
##'
##' @name rpf.sample
##' @param theta either a vector (for 1 dimension) or a matrix (for >1
##' dimension) of person abilities or the number of response patterns
##' to generate randomly
##' @param items a list of item models
##' @param params a list or matrix of item parameters. If omitted, random item
##' parameters are generated for each item model.
##' @param design a matrix assigning person abilities to item factors
##' @param prefix Column names are taken from param or items.
##' If no column names are available, some will be generated using
##' the given prefix.
##' @param mean mean vector of latent distribution (optional)
##' @param cov covariance matrix of latent distribution (optional)
##' @return Returns a data frame of response patterns
##' @export
##' @examples
##' # 1 dimensional items
##' i1 <- rpf.drm()
##' i1.p <- rpf.rparam(i1)
##' i2 <- rpf.nrm(outcomes=3)
##' i2.p <- rpf.rparam(i2)
##' rpf.sample(5, list(i1,i2), list(i1.p, i2.p))
##'
##' # multidimensional items
##' numItems <- 4
##' items <- vector("list", numItems)
##' correct <- vector("list", numItems)
##'
##' i1 <- rpf.drm(factors=2)
##' i2 <- rpf.drm(factors=1, multidimensional=TRUE)
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
##' @references
##' Cai, L. (2010). A two-tier full-information item factor analysis
##' model with applications. \emph{Psychometrika, 75}, 581-612.
rpf.sample <- function(theta, items, params, design, prefix="i",
                       mean=NULL, cov=NULL) {
  numItems <- length(items)
  maxDim <- max(vapply(items, function(i) i@factors, 0))
  if (missing(design)) {
    if (maxDim > 1) {
      design <- matrix(rep(1:maxDim, numItems), nrow=maxDim)
      design[sapply(items, function(i) 1:maxDim > i@factors)] <- NA
    } else {
      design <- matrix(rep(1, numItems), nrow=1)
    }
  }
  maxAbilities <- max(design, na.rm=TRUE)

  if (maxDim > 1 && any(dim(design) != c(maxDim,numItems))) {
    stop(paste("The design matrix must have", maxDim, "rows and ",numItems,"columns"))
  }

  numPeople <- NA
  if (is.numeric(theta) && length(theta) == 1) {
    if (theta <= 1) stop("Request at least 2 samples")
    numPeople <- theta
    if (missing(mean)) mean <- rep(0, maxAbilities)
    if (length(mean) != maxAbilities) stop(paste("Mean vector must have length",maxAbilities))
    if (missing(cov)) cov <- diag(maxAbilities)
    if (any(dim(cov) != maxAbilities)) stop(paste("Cov matrix must be square matrices of size",maxAbilities))
    theta <- array(t(rmvnorm(numPeople, mean=mean, sigma=cov)),
                   dim=c(maxAbilities, numPeople))
  } else if (maxDim == 1 && is.vector(theta)) {
    numPeople <- length(theta)
    theta <- array(theta, dim=c(maxAbilities, numPeople))
  } else {
      if (dim(theta)[1] > maxAbilities) {
          stop(paste("Only", maxAbilities, "abilities but theta provides", dim(theta)[1],
                     "-- maybe transpose theta?"))
      }
    numPeople <- dim(theta)[2]
  }

  if (missing(params)) {
    params <- lapply(items, rpf.rparam)
  }

  outcomes <- vapply(items, function(i) i@outcomes, 0)
  
  ret <- list()
  for (ix in 1:numItems) {
    i <- items[[ix]]
    param <- c()
    if (is.list(params)) {
      param <- params[[ix]]
    } else {
      param <- params[,ix]  # item parameters are in columns
    }

    if (length(param) < rpf.numParam(i)) {
      stop(paste("Item",class(i),"needs",rpf.numParam(i),
                 "parameters but only",length(param),"given"))
    }

    cols <- design[,ix]
    cols <- cols[!is.na(cols)]
    i.theta <- theta[cols,,drop=FALSE]
    P <- rpf.prob(i, param[1:rpf.numParam(i)], i.theta)
#    if (any(is.na(P))) stop(paste("Item", i@spec, "with param", param," produced NAs"))
    ret1 <- apply(P, 2, sample, x=1:i@outcomes, size=1, replace=F)
    ret1 <- factor(ret1, levels=1:i@outcomes, ordered=TRUE)
    attr(ret1, 'mxFactor') <- TRUE  # for OpenMx
    ret[[ix]] <- ret1
  }
  ret <- as.data.frame(ret)
  name <- colnames(params)
  if (is.null(name)) name <- names(items)
  if (is.null(name)) name <- paste(prefix,1:numItems,sep="")
  colnames(ret) <- name
  return(ret)
}

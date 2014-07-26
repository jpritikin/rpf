##' Convert an OpenMx MxModel object into an IFA group
##'
##' @param mxModel MxModel object
##' @param data observed data (otherwise the data will be taken from the mxModel)
##' @return a groups with item parameters and latent distribution
as.IFAgroup <- function(mxModel, data=NULL) {
  if (!is(mxModel$expectation, "MxExpectationBA81")) {
    stop(paste("Don't know how to create an IFA group from",
               class(mxModel$expectation)))
  }

  mat <- mxModel$expectation$item
  if (length(grep("\\.", mat))) {
    stop(paste("Don't know how to obtain the item matrix", mat))
  }
  itemMat <- mxModel[[mat]]
  if (is.null(itemMat)) {
    stop(paste("Item matrix", mat, "not found"))
  }
  
  ret <- list(spec = mxModel$expectation$ItemSpec,
              param = itemMat$values,
              free = itemMat$free,
              qpoints = mxModel$expectation$qpoints,
              qwidth = mxModel$expectation$qwidth)
  
  mat <- mxModel$expectation$mean
  if (length(grep("\\.", mat))) {
    stop(paste("Don't know how to obtain the mean matrix", mat))
  }
  meanMat <- mxModel[[mat]]
  if (!is.null(meanMat)) {
    ret$mean <- meanMat$values
    ret$meanFree <- meanMat$free
  }

  mat <- mxModel$expectation$cov
  if (length(grep("\\.", mat))) {
    stop(paste("Don't know how to obtain the cov matrix", mat))
  }
  covMat <- mxModel[[mat]]
  if (!is.null(covMat)) {
    ret$cov <- covMat$values
    ret$covFree <- covMat$free
  }
  
  if (!missing(data)) {
    ret$data <- data
  } else if (!is.null(mxModel$data)) {
    mxData <- mxModel$data
    if (mxData$type != "raw") {
      stop(paste("Not sure how to handle data of type", mxData$type))
    }
    if (mxData$.isSorted) {
      unsort <- match(0:(length(mxData$indexVector)-1), mxData$indexVector)
      ret$data <- mxData$observed[unsort,]
    } else {
      ret$data <- mxData$observed
    }
  }

  if (!is.na(mxModel$expectation$minItemsPerScore)) {
    ret$minItemsPerScore <- mxModel$expectation$minItemsPerScore
  }
  if (!is.na(mxModel$expectation$weightColumn)) {
    ret$weightColumn <- mxModel$expectation$weightColumn
  }
  ret
}

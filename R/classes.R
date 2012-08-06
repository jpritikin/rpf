##' rpf - Response Probability Functions
##'
##' The purpose of this package is to factor out logic and math common
##' to Item Response Theory fitting, diagnostics, and analysis.  It is
##' envisioned as core support code suitable for more specialized IRT
##' packages to build upon.
##'
##' This package provides optimized, low-level functions to map
##' parameters to response probabilities for dichotomous (1PL, 2PL and
##' 3PL) \code{\link{rpf.drm}} and polytomous (graded response,
##' partial credit/generalized partial credit \code{\link{rpf.gpcm}},
##' nominal, and multiple-choice model) items. Both unidimensional
##' and multidimensional versions of the models will be available.
##'
##' This package could also accrete functions to support plotting (but
##' not the actual plot functions).
##'
##' @section Warning:
##' The API is not stable at this time. You have been warned.
##' 
##' @docType package
##' @rdname rpf.introduction
##' @name An introduction
##' @seealso TODO in the source distribution
NULL

##' The base class for response probability functions.
##' @name Class rpf.base
##' @rdname rpf.base-class
##' @aliases rpf.base-class
##' @export
setClass("rpf.base",
         representation(numOutcomes="numeric", "VIRTUAL"))

##' Map an item model, item parameters, and person trait score into a
##' probability vector
##'
##' @param m an item model
##' @param param item parameters
##' @param theta a single person's trait score
##' @return a vector of probabilities
##' @docType methods
##' @aliases
##' rpf.prob,rpf.drm,numeric,numeric-method
##' rpf.prob,rpf.gpcm,numeric,numeric-method
##' @export
##' @seealso
##' See \code{\link{rpf.drm}} and \code{\link{rpf.gpcm}} to create item models.
##' 
##' See \code{\link{rpf.rparam}} to create item parameters.
##' @examples
##' i1 <- rpf.drm()
##' i1.p <- rpf.rparam(i1)
##' rpf.prob(i1, c(i1.p), -1)   # low trait score
##' rpf.prob(i1, c(i1.p), 0)    # average trait score
##' rpf.prob(i1, c(i1.p), 1)    # high trait score
setGeneric("rpf.prob", function(m, param, theta) standardGeneric("rpf.prob"))

##' Log likelihood of item parameters with respect to Bayesian prior
##'
##' This function calculates the log likelihood of the item parameters
##' with respect to the item model's Bayesian prior. This function is
##' typically used by model fitting algorithms.
##' 
##' @param m an item model
##' @param param item parameters
##' @docType methods
##' @aliases
##' rpf.logLik,rpf.drm,numeric-method
##' rpf.logLik,rpf.gpcm,numeric-method
##' @export
##' @return a log likelihood (not -2 * log likelihood)
##' @examples
##' i1 <- rpf.drm()
##' i1.p <- rpf.rparam(i1)
##' rpf.logLik(i1, c(i1.p))
setGeneric("rpf.logLik", function(m, param) standardGeneric("rpf.logLik"))

##' Dimensions of an item's parameters
##'
##' This method is available in case you want to dimension a matrix to
##' exactly the right size to hold item parameters.
##' 
##' @param m an item model
##' @return a vector of dimension sizes
##' @docType methods
##' @aliases
##' rpf.paramDim,rpf.drm-method
##' rpf.paramDim,rpf.gpcm-method
##' @export
##' @examples
##' i1 <- rpf.drm()
##' i2 <- rpf.gpcm(numOutcomes=5)
##' apply(sapply(list(i1,i2), rpf.paramDim), 1, max)
setGeneric("rpf.paramDim", function(m) standardGeneric("rpf.paramDim"))

##' Generates item parameters
##'
##' This function generates random item parameters taking the Bayesian
##' priors into account.
##' 
##' @param m an item model
##' @return item parameters
##' @docType methods
##' @aliases
##' rpf.rparam,rpf.drm-method
##' rpf.rparam,rpf.gpcm-method
##' @export
##' @examples
##' i1 <- rpf.drm()
##' rpf.rparam(i1)
setGeneric("rpf.rparam", function(m) standardGeneric("rpf.rparam"))

##' Initial item parameters
##'
##' This function generates item parameters suitable for use as
##' initial starting values in an IRT fitting algorithm.
##' 
##' @param m an item model
##' @return item parameters
##' @docType methods
##' @aliases
##' rpf.startingParam,rpf.drm-method
##' rpf.startingParam,rpf.gpcm-method
##' @export
##' @examples
##' i1 <- rpf.drm()
##' rpf.startingParam(i1)
setGeneric("rpf.startingParam", function(m) standardGeneric("rpf.startingParam"))

##' Extract location related item parameters
##'
##' @param m an item model
##' @param param item parameters
##' @return location related item parameters
##' @docType methods
##' @aliases
##' rpf.getLocation,rpf.drm,numeric-method
##' rpf.getLocation,rpf.gpcm,numeric-method
##' @export
##' @seealso \code{\link{rpf.setLocation}}
setGeneric("rpf.getLocation", function(m,param) standardGeneric("rpf.getLocation"))

##' Set location related item parameters
##'
##' @param m an item model
##' @param param item parameters
##' @param loc location parameters, typically from \code{\link{rpf.getLocation}}
##' @return updated item parameters
##' @docType methods
##' @aliases
##' rpf.setLocation,rpf.drm,numeric,numeric-method
##' rpf.setLocation,rpf.gpcm,numeric,numeric-method
##' @export
##' @seealso \code{\link{rpf.getLocation}}
setGeneric("rpf.setLocation", function(m,param,loc) standardGeneric("rpf.setLocation"))

##' The ogive constant
##'
##' Models built on the logistic function take an argument \code{D}
##' where you can pass in the ogive constant to obtain a response
##' curve very similar to the Normal cumulative distribution function.
##' In recent years, the logistic has grown in favor, and therefore,
##' \code{D} defaults to 1 in this package (Baker & Kim, 2004, pp. 14-18).
##' 
##' @export
##' @references Baker & Kim (2004). Item Response Theory: Parameter
##' Estimation Techniques. Marcel Dekker, Inc.
##' 
rpf.ogive <- 1.702

##' The base class for 1 dimensional response probability functions.
##' @name Class rpf.1dim
##' @rdname rpf.1dim-class
##' @aliases rpf.1dim-class
##' @export
setClass("rpf.1dim", contains='rpf.base',
         representation("VIRTUAL"))

##' The base class for multi-dimensional response probability functions.
##' @name Class rpf.mdim
##' @rdname rpf.mdim-class
##' @aliases rpf.mdim-class
##' @export
setClass("rpf.mdim", contains='rpf.base',
         representation(dimensions="numeric",
                        "VIRTUAL"))

##' The unidimensional generalized partial credit item model.
##'
##' @export
##' @name Class rpf.gpcm
##' @rdname rpf.gpcm-class
##' @aliases rpf.gpcm-class
setClass("rpf.gpcm", contains='rpf.1dim',
         representation(D="numeric",
                        a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric"))

##' Unidimensional dichotomous item models (1PL, 2PL, and 3PL).
##'
##' @export
##' @name Class rpf.drm
##' @rdname rpf.drm-class
##' @aliases rpf.drm-class
setClass("rpf.drm", contains='rpf.1dim',
         representation(D="numeric",
                        guessing="numeric",
                        a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric",
                        c.prior.alpha="numeric",
                        c.prior.beta="numeric"))

##' Randomly sample response patterns given a list of items
##'
##' Returns a random sample of response patterns given
##' a list of item models and parameters.
##'
##' @name rpf.sample
##' @usage
##' rpf.sample(theta, items, params)
##' @param theta either a vector of trait abilities or
##' the number of abilities to draw from N(0,1)
##' @param items a list of item models
##' @param params a list of item parameters. If omitted, random item
##' parameters are generated for each item model.
##' @return Returns a matrix of response patterns
##' @export
##' @examples
##' i1 <- rpf.drm()
##' i1.p <- rpf.rparam(i1)
##' i2 <- rpf.gpcm(numOutcomes=3)
##' i2.p <- rpf.rparam(i2)
##' rpf.sample(5, list(i1,i2), list(i1.p, i2.p))
##' @seealso \code{\link{sample}}
rpf.sample <- function(theta, items, params) {
    if (length(theta) == 1) {
        theta <- rnorm(theta)
    }
    numPeople <- length(theta)
    numItems <- length(items)
    if (missing(params)) {
        params <- lapply(items, rpf.rparam)
    }
    outcomes <- vapply(items, function(i) i@numOutcomes, 0)
    
    ret <- array(dim=c(numPeople, numItems))
    for (ix in 1:numItems) {
        i <- items[[ix]]
        param <- params[[ix]]
        P <- rpf.prob(i, c(param), theta)
        ret[,ix] <- apply(P, c(1), sample, x=1:i@numOutcomes, size=1, replace=F)
    }
    return(ret)
}

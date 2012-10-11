##' rpf - Response Probability Functions
##'
##' The purpose of this package is to factor out logic and math common
##' to Item Response Theory fitting, diagnostics, and analysis.  It is
##' envisioned as core support code suitable for more specialized IRT
##' packages to build upon.
##'
##' This package provides optimized, low-level functions to map
##' parameters to response probabilities for dichotomous (1PL, 2PL and
##' 3PL) \code{\link{rpf.drm}} and polytomous (graded response
##' \code{\link{rpf.grm}}, partial credit/generalized partial credit
##' \code{\link{rpf.gpcm}}, nominal \code{\link{rpf.nrm}}, and
##' multiple-choice model \code{\link{rpf.mcm}}) items. Both
##' unidimensional and multidimensional versions of the models are
##' available.
##'
##' Item model parameters are passed around as a numeric vector. A 1D
##' matrix is also acceptable. Regardless of model, parameters are
##' always ordered as follows: discrimination ("a"), difficulty ("b"),
##' and guessing ("c").
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
##' @useDynLib rpf
##' @seealso
##' See \code{\link{rpf.rparam}} and \code{\link{rpf.startingParam}}
##' to create item parameters.
NULL

##' The base class for response probability functions.
##' @name Class rpf.base
##' @rdname rpf.base-class
##' @aliases rpf.base-class
##' @export
setClass("rpf.base",
         representation(numOutcomes="numeric",
                        numParam="numeric",
                        dimensions="numeric",
                        "VIRTUAL"))

##' Map an item model, item parameters, and person trait score into a
##' probability vector
##'
##' @param m an item model
##' @param param item parameters
##' @param theta the trait score(s)
##' @return a vector of probabilities. For dichotomous items,
##' probabilities are returned in the order incorrect, correct.
##' Although redundent, both incorrect and correct probabilities are
##' returned for API consistency with polytomous item models.
##' @docType methods
##' @aliases
##' rpf.prob,rpf.1dim.drm,numeric,numeric-method
##' rpf.prob,rpf.mdim.drm,numeric,matrix-method
##' rpf.prob,rpf.1dim.grm,numeric,numeric-method
##' rpf.prob,rpf.mdim.grm,numeric,numeric-method
##' rpf.prob,rpf.1dim.gpcm,numeric,numeric-method
##' rpf.prob,rpf.mdim.gpcm,numeric,matrix-method
##' rpf.prob,rpf.mdim.nrm,numeric,matrix-method
##' rpf.prob,rpf.mdim.mcm,numeric,matrix-method
##' rpf.prob,rpf.mdim.grm,numeric,matrix-method
##' @export
##' @examples
##' i1 <- rpf.drm()
##' i1.p <- rpf.rparam(i1)
##' rpf.prob(i1, c(i1.p), -1)   # low trait score
##' rpf.prob(i1, c(i1.p), c(0,1))    # average and high trait score
setGeneric("rpf.prob", function(m, param, theta) standardGeneric("rpf.prob"))

##' Turn a matrix of parameters into a vector for
##' \code{\link{rpf.prob}}.
##' 
##' @name rpf.prob wrapper1
##' @rdname rpf.prob.wrapper1
##' @aliases rpf.prob,rpf.base,matrix,numeric-method
##' @docType methods
setMethod("rpf.prob", signature(m="rpf.base", param="matrix",
                                theta="numeric"),
          function(m, param, theta) {
            rpf.prob(m, as.numeric(param), theta)
          })

##' Turn a matrix of parameters into a vector for
##' \code{\link{rpf.prob}}.
##' 
##' @name rpf.prob wrapper2
##' @rdname rpf.prob.wrapper2
##' @aliases rpf.prob,rpf.base,matrix,matrix-method
##' @docType methods
setMethod("rpf.prob", signature(m="rpf.base", param="matrix",
                                theta="matrix"),
          function(m, param, theta) {
            rpf.prob(m, as.numeric(param), theta)
          })

##' Map an item model, item parameters, and person trait score into a
##' information vector
##'
##' @param m an item model
##' @param param item parameters
##' @param theta the trait score(s)
##' @return Fisher information
##' @docType methods
##' @aliases
##' rpf.info,rpf.1dim.drm,numeric,numeric-method
##' @export
##' @examples
##' i1 <- rpf.drm()
##' i1.p <- c(.6,1,.1)
##' theta <- seq(0,3,.05)
##' plot(theta, rpf.info(i1, i1.p, theta), type="l")
setGeneric("rpf.info", function(m, param, theta) standardGeneric("rpf.info"))

##' Generates item parameters
##'
##' This function generates random item parameters taking the Bayesian
##' priors into account.
##' 
##' @param m an item model
##' @return item parameters
##' @docType methods
##' @aliases
##' rpf.rparam,rpf.1dim.drm-method
##' rpf.rparam,rpf.mdim.drm-method
##' rpf.rparam,rpf.1dim.graded-method
##' rpf.rparam,rpf.mdim.graded-method
##' rpf.rparam,rpf.mdim.nrm-method
##' rpf.rparam,rpf.mdim.mcm-method
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
##' rpf.startingParam,rpf.1dim.drm-method
##' rpf.startingParam,rpf.mdim.drm-method
##' rpf.startingParam,rpf.1dim.graded-method
##' rpf.startingParam,rpf.mdim.graded-method
##' rpf.startingParam,rpf.mdim.nrm-method
##' rpf.startingParam,rpf.mdim.mcm-method
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
##' rpf.getLocation,rpf.1dim.drm,numeric-method
##' rpf.getLocation,rpf.mdim.drm,numeric-method
##' rpf.getLocation,rpf.1dim.graded,numeric-method
##' rpf.getLocation,rpf.mdim.graded,numeric-method
##' rpf.getLocation,rpf.mdim.nrm,numeric-method
##' rpf.getLocation,rpf.mdim.mcm,numeric-method
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
##' rpf.setLocation,rpf.1dim.drm,numeric,numeric-method
##' rpf.setLocation,rpf.mdim.drm,numeric,numeric-method
##' rpf.setLocation,rpf.1dim.graded,numeric,numeric-method
##' rpf.setLocation,rpf.mdim.graded,numeric,numeric-method
##' rpf.setLocation,rpf.mdim.nrm,numeric,numeric-method
##' rpf.setLocation,rpf.mdim.mcm,numeric,numeric-method
##' @export
##' @seealso \code{\link{rpf.getLocation}}
setGeneric("rpf.setLocation", function(m,param,loc) standardGeneric("rpf.setLocation"))

##' The ogive constant
##'
##' Models built on the logistic function take an argument \code{D}
##' where you can pass in the ogive constant to obtain a response
##' curve very similar to the Normal cumulative distribution function
##' (Haley, 1952).
##' In recent years, the logistic has grown in favor, and therefore,
##' \code{D} defaults to 1 in this package (Baker & Kim, 2004, pp. 14-18).
##' 
##' @export
##' @references Baker & Kim (2004). Item Response Theory: Parameter
##' Estimation Techniques. Marcel Dekker, Inc.
##'
##' Haley, D. C. (1952). Estimation of the dosage mortality
##' relationship when the dose is subject to error (Technical Report
##' No. 15). Stanford University Applied Mathematics and Statistics
##' Laboratory, Stanford, CA.
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
         representation("VIRTUAL"))

##' The base class for 1 dimensional graded response probability functions.
##' 
##' This class contains methods common to both the generalized partial
##' credit model and the graded response model.
##'
##' @name Class rpf.1dim.graded
##' @rdname rpf.1dim.graded-class
##' @aliases rpf.1dim.graded-class
##' @export
setClass("rpf.1dim.graded", contains='rpf.1dim',
         representation("VIRTUAL"))

##' The base class for multi-dimensional graded response probability
##' functions.
##'
##' This class contains methods common to both the generalized partial
##' credit model and the graded response model.
##'
##' @name Class rpf.mdim.graded
##' @rdname rpf.mdim.graded-class
##' @aliases rpf.mdim.graded-class
##' @export
setClass("rpf.mdim.graded", contains='rpf.mdim',
         representation("VIRTUAL"))

##' The unidimensional graded response item model.
##'
##' @export
##' @name Class rpf.1dim.grm
##' @rdname rpf.1dim.grm-class
##' @aliases rpf.1dim.grm-class
setClass("rpf.1dim.grm", contains='rpf.1dim.graded',
         representation(D="numeric",
                        a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric"))

##' The unidimensional generalized partial credit item model.
##'
##' @export
##' @name Class rpf.1dim.gpcm
##' @rdname rpf.1dim.gpcm-class
##' @aliases rpf.1dim.gpcm-class
setClass("rpf.1dim.gpcm", contains='rpf.1dim.graded',
         representation(D="numeric",
                        a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric"))

##' Unidimensional dichotomous item models (1PL, 2PL, and 3PL).
##'
##' @export
##' @name Class rpf.1dim.drm
##' @rdname rpf.1dim.drm-class
##' @aliases rpf.1dim.drm-class
setClass("rpf.1dim.drm", contains='rpf.1dim',
         representation(D="numeric",
                        guessing="numeric",
                        a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric",
                        c.prior.alpha="numeric",
                        c.prior.beta="numeric"))

##' Multidimensional dichotomous item models (M1PL, M2PL, and M3PL).
##'
##' @export
##' @name Class rpf.mdim.drm
##' @rdname rpf.mdim.drm-class
##' @aliases rpf.mdim.drm-class
setClass("rpf.mdim.drm", contains='rpf.mdim',
         representation(D="numeric",
                        guessing="numeric",
                        a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric",
                        c.prior.alpha="numeric",
                        c.prior.beta="numeric"))

##' The multidimensional graded response item model.
##'
##' @export
##' @name Class rpf.mdim.grm
##' @rdname rpf.mdim.grm-class
##' @aliases rpf.mdim.grm-class
setClass("rpf.mdim.grm", contains='rpf.mdim.graded',
         representation(D="numeric",
                        a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric"))

##' The multidimensional generalized partial credit item model.
##'
##' @export
##' @name Class rpf.mdim.gpcm
##' @rdname rpf.mdim.gpcm-class
##' @aliases rpf.mdim.gpcm-class
setClass("rpf.mdim.gpcm", contains='rpf.mdim.graded',
         representation(D="numeric",
                        a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric"))

##' The nominal response item model (both unidimensional and
##' multidimensional models have the same parameterization).
##'
##' @export
##' @name Class rpf.mdim.nrm
##' @rdname rpf.mdim.nrm-class
##' @aliases rpf.mdim.nrm-class
setClass("rpf.mdim.nrm", contains='rpf.mdim',
         representation(a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric"))

##' The multiple-choice response item model (both unidimensional and
##' multidimensional models have the same parameterization).
##'
##' @export
##' @name Class rpf.mdim.mcm
##' @rdname rpf.mdim.mcm-class
##' @aliases rpf.mdim.mcm-class
setClass("rpf.mdim.mcm", contains='rpf.mdim',
         representation(a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric",
                        c.prior.alpha="numeric",
                        c.prior.beta="numeric"))

##' Randomly sample response patterns given a list of items
##'
##' Returns a random sample of response patterns given
##' a list of unidimensional item models and parameters.
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
  if (max(vapply(items, function(i) i@dimensions, 0)) > 1) {
    stop("Can only sample unidimensional items")
  }
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
    P <- rpf.prob(i, param, theta)
    ret[,ix] <- apply(P, c(1), sample, x=1:i@numOutcomes, size=1, replace=F)
  }
  return(ret)
}

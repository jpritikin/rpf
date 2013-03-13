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
##' and guessing ("c"). If person ability ranges from low negative to
##' high positive then probabilities are output from incorrect to
##' correct. That is, a low ability person (e.g., ability = -2) will
##' be more likely to get an item incorrect than correct. For example,
##' a dichotomous model that returns [.25, .75] indicates a
##' probability of .25 for incorrect and .75 for correct.  A
##' polytomous model will have the most incorrect probability at index
##' 1 and the most correct probability at the maximum index.
##' 
##' All models are always in the logistic metric. To obtain normal
##' ogive discrimination parameters, divide slope parameters by
##' \code{\link{rpf.ogive}}. Item models are estimated in
##' slope-intercept form unless the traditional parameterization is
##' specifically requested.
##'
##' This package could also accrete functions to support plotting (but
##' not the actual plot functions).
##'
##' @section Warning: The API is not stable at this time. You have
##' been warned. In particular, I anticipating transposing some of the
##' input and output matrices.
##' 
##' @docType package
##' @rdname rpf.introduction
##' @name An introduction
##' @useDynLib rpf
##' @references Thissen, D. and Steinberg, L. (1986). A taxonomy of
##' item response models. \emph{Psychometrika 51}(4), 567-577.
##' @seealso
##' See \code{\link{rpf.rparam}} to create item parameters.
NULL

##' The base class for response probability functions.
##'
##' Item specifications should not be modified after creation.
##' 
##' @name Class rpf.base
##' @rdname rpf.base-class
##' @aliases rpf.base-class
##' @export
setClass("rpf.base",
         representation(spec="numeric",
                        numOutcomes="numeric",
                        numParam="numeric",
                        dimensions="numeric",
                        "VIRTUAL"))

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

##' Length of the item model vector
##' @aliases
##' rpf.numSpec,rpf.base-method
##' rpf_numSpec_wrapper
setGeneric("rpf.numSpec", function(m) standardGeneric("rpf.numSpec"))

setMethod("rpf.numSpec", signature(m="rpf.base"),
          function(m) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              .Call(rpf_numSpec_wrapper, m@spec)
            }
          })

##' Length of the item parameter vector
##' @aliases
##' rpf.numParam,rpf.base-method
##' rpf_numParam_wrapper
setGeneric("rpf.numParam", function(m) standardGeneric("rpf.numParam"))

setMethod("rpf.numParam", signature(m="rpf.base"),
          function(m) {
            if (length(m@spec)==0) {
              m@numParam
            } else {
              .Call(rpf_numParam_wrapper, m@spec)
            }
          })

##' Log likelihood of the item model parameters given the Bayesian prior
##' @aliases
##' rpf.prior,rpf.base,numeric-method
##' rpf_prior_wrapper
setGeneric("rpf.prior", function(m, param) standardGeneric("rpf.prior"))

setMethod("rpf.prior", signature(m="rpf.base", param="numeric"),
          function(m, param) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              .Call(rpf_prior_wrapper, m@spec, param)
            }
          })

##' Item parameter gradients
##'
##' Evaluate the partial derivatives of the log likelihood with
##' respect to each parameter at \code{where} with \code{weight}.
##'
##' @param m item model
##' @param param item parameters
##' @param where location in the latent space
##' @param weight per outcome weights (typically derived by observation)
##' @return derivative of the log likelihood with respect to each parameter evaluated at \code{where}
##' @aliases
##' rpf.gradient,rpf.base,numeric,numeric,numeric-method
##' rpf_gradient_wrapper
setGeneric("rpf.gradient", function(m, param, where, weight) standardGeneric("rpf.gradient"))

setMethod("rpf.gradient", signature(m="rpf.base", param="numeric",
                                    where="numeric", weight="numeric"),
          function(m, param, where, weight) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              .Call(rpf_gradient_wrapper, m@spec, param, where, weight)
            }
          })

##' Map an item model, item parameters, and person trait score into a
##' probability vector
##'
##' In some cases, this function is implemented in terms of \code{\link{rpf.logprob}}.
##'
##' @param m an item model
##' @param param item parameters
##' @param theta the trait score(s)
##' @return a vector of probabilities. For dichotomous items,
##' probabilities are returned in the order incorrect, correct.
##' Although redundent, both incorrect and correct probabilities are
##' returned in the dichotomous case for API consistency with
##' polytomous item models.
##' @docType methods
##' @aliases
##' rpf.prob,rpf.1dim,numeric,numeric-method
##' rpf.prob,rpf.mdim,numeric,numeric-method
##' rpf.prob,rpf.mdim,numeric,matrix-method
##' rpf.prob,rpf.base,data.frame,numeric-method
##' rpf.prob,rpf.base,matrix,numeric-method
##' rpf.prob,rpf.base,matrix,matrix-method
##' rpf.prob,rpf.1dim,numeric,matrix-method
##' rpf.prob,rpf.1dim.grm,numeric,numeric-method
##' rpf.prob,rpf.mdim.grm,numeric,numeric-method
##' rpf.prob,rpf.mdim.gpcm,numeric,matrix-method
##' rpf.prob,rpf.mdim.nrm,numeric,matrix-method
##' rpf.prob,rpf.mdim.mcm,numeric,matrix-method
##' rpf.prob,rpf.mdim.grm,numeric,matrix-method
##' rpf_prob_wrapper
##' @export
##' @examples
##' i1 <- rpf.drm()
##' i1.p <- rpf.rparam(i1)
##' rpf.prob(i1, c(i1.p), -1)   # low trait score
##' rpf.prob(i1, c(i1.p), c(0,1))    # average and high trait score
setGeneric("rpf.prob", function(m, param, theta) standardGeneric("rpf.prob"))

##' Map an item model, item parameters, and person trait score into a
##' probability vector
##'
##' @param m an item model
##' @param param item parameters
##' @param theta the trait score(s)
##' @return a vector of probabilities. For dichotomous items,
##' probabilities are returned in the order incorrect, correct.
##' Although redundent, both incorrect and correct probabilities are
##' returned in the dichotomous case for API consistency with
##' polytomous item models.
##' @docType methods
##' @aliases
##' rpf.logprob,rpf.1dim,numeric,numeric-method
##' rpf.logprob,rpf.1dim,numeric,matrix-method
##' rpf.logprob,rpf.mdim,numeric,matrix-method
##' rpf.logprob,rpf.mdim,numeric,numeric-method
##' rpf_logprob_wrapper
##' @export
##' @examples
##' i1 <- rpf.gpcm()
##' i1.p <- rpf.rparam(i1)
##' rpf.logprob(i1, c(i1.p), -1)   # low trait score
##' rpf.logprob(i1, c(i1.p), c(0,1))    # average and high trait score
setGeneric("rpf.logprob", function(m, param, theta) standardGeneric("rpf.logprob"))

setMethod("rpf.logprob", signature(m="rpf.1dim", param="numeric", theta="numeric"),
          function(m, param, theta) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              t(.Call(rpf_logprob_wrapper, m@spec, param, theta))
            }
          })

setMethod("rpf.logprob", signature(m="rpf.mdim", param="numeric", theta="matrix"),
          function(m, param, theta) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              t(.Call(rpf_logprob_wrapper, m@spec, param, t(theta)))
            }
          })

setMethod("rpf.logprob", signature(m="rpf.mdim", param="numeric", theta="numeric"),
          function(m, param, theta) {
            rpf.logprob(m, param, as.matrix(theta))
          })

setMethod("rpf.logprob", signature(m="rpf.1dim", param="numeric", theta="matrix"),
          function(m, param, theta) {
            rpf.logprob(m, param, as.numeric(theta))
          })

setMethod("rpf.prob", signature(m="rpf.1dim", param="numeric", theta="numeric"),
          function(m, param, theta) {
            if (length(m@spec)==0) {
              exp(rpf.logprob(m, param, theta))
            } else {
              t(.Call(rpf_prob_wrapper, m@spec, param, theta))
            }
          })

setMethod("rpf.prob", signature(m="rpf.mdim", param="numeric", theta="matrix"),
          function(m, param, theta) {
            if (length(m@spec)==0) {
              exp(rpf.logprob(m, param, theta))
            } else {
              t(.Call(rpf_prob_wrapper, m@spec, param, t(theta)))
            }
          })

setMethod("rpf.prob", signature(m="rpf.base", param="data.frame", theta="numeric"),
          function(m, param, theta) {
            exp(rpf.logprob(m, as.numeric(param), theta))
          })

setMethod("rpf.prob", signature(m="rpf.base", param="matrix", theta="numeric"),
          function(m, param, theta) {
            rpf.prob(m, as.numeric(param), theta)
          })

setMethod("rpf.prob", signature(m="rpf.base", param="matrix", theta="matrix"),
          function(m, param, theta) {
            rpf.prob(m, as.numeric(param), theta)
          })

setMethod("rpf.prob", signature(m="rpf.1dim", param="numeric", theta="matrix"),
          function(m, param, theta) {
            rpf.prob(m, param, as.numeric(theta))
          })

setMethod("rpf.prob", signature(m="rpf.mdim", param="numeric", theta="numeric"),
          function(m, param, theta) {
            rpf.prob(m, param, as.matrix(theta))
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
##' rpf.info,rpf.base,data.frame,numeric-method
##' rpf.info,rpf.1dim.drm,numeric,numeric-method
##' rpf.info,rpf.mdim.drm,numeric,numeric-method
##' rpf.info,rpf.1dim.graded,numeric,numeric-method
##' @export
##' @examples
##' i1 <- rpf.drm()
##' i1.p <- c(.6,1,.1)
##' theta <- seq(0,3,.05)
##' plot(theta, rpf.info(i1, i1.p, theta), type="l")
##' @references
##' Muraki, E. (1993) Information functions of the generalized partial credit model.
##' Applied Psychological Measurement, 17(4), 351-363.
setGeneric("rpf.info", function(m, param, theta) standardGeneric("rpf.info"))

setMethod("rpf.info", signature(m="rpf.base", param="data.frame", theta="numeric"),
          function(m, param, theta) {
            rpf.info(m, as.numeric(param), theta)
          })

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

##' The ogive constant
##'
##' The ogive constant can be multiplied by the discrimination
##' parameter to obtain a response curve very similar to the Normal
##' cumulative distribution function (Haley, 1952; Molenaar, 1974).
##' In recent years, the logistic has grown in favor, and therefore,
##' this package does not offer any special support for this
##' transformation (Baker & Kim, 2004, pp. 14-18).
##' 
##' @export
##' @references
##' Baker & Kim (2004). \emph{Item Response Theory: Parameter
##' Estimation Techniques.} Marcel Dekker, Inc.
##'
##' Haley, D. C. (1952). \emph{Estimation of the dosage mortality
##' relationship when the dose is subject to error} (Technical Report
##' No. 15). Stanford University Applied Mathematics and Statistics
##' Laboratory, Stanford, CA.
##'
##' Molenaar, W. (1974). De logistische en de normale kromme [The
##' logistic and the normal curve]. \emph{Nederlands Tijdschrift voor de
##' Psychologie} 29, 415-420.
rpf.ogive <- 1.702

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
         representation(a.prior.sdlog="numeric"))

##' The unidimensional generalized partial credit item model.
##'
##' @export
##' @name Class rpf.1dim.gpcm
##' @rdname rpf.1dim.gpcm-class
##' @aliases rpf.1dim.gpcm-class
setClass("rpf.1dim.gpcm", contains='rpf.1dim.graded')

##' Unidimensional dichotomous item models (1PL, 2PL, and 3PL).
##'
##' @export
##' @name Class rpf.1dim.drm
##' @rdname rpf.1dim.drm-class
##' @aliases rpf.1dim.drm-class
setClass("rpf.1dim.drm", contains='rpf.1dim',
         representation(c.prior.logit="numeric"))

##' Multidimensional dichotomous item models (M1PL, M2PL, and M3PL).
##'
##' @export
##' @name Class rpf.mdim.drm
##' @rdname rpf.mdim.drm-class
##' @aliases rpf.mdim.drm-class
setClass("rpf.mdim.drm", contains='rpf.mdim',
         representation(c.prior.logit="numeric"))

##' The multidimensional graded response item model.
##'
##' @export
##' @name Class rpf.mdim.grm
##' @rdname rpf.mdim.grm-class
##' @aliases rpf.mdim.grm-class
setClass("rpf.mdim.grm", contains='rpf.mdim.graded',
         representation(a.prior.sdlog="numeric"))

##' The multidimensional generalized partial credit item model.
##'
##' @export
##' @name Class rpf.mdim.gpcm
##' @rdname rpf.mdim.gpcm-class
##' @aliases rpf.mdim.gpcm-class
setClass("rpf.mdim.gpcm", contains='rpf.mdim.graded')

##' The nominal response item model (both unidimensional and
##' multidimensional models have the same parameterization).
##'
##' @export
##' @name Class rpf.mdim.nrm
##' @rdname rpf.mdim.nrm-class
##' @aliases rpf.mdim.nrm-class
setClass("rpf.mdim.nrm", contains='rpf.mdim',
         representation(a.prior.sdlog="numeric"))

##' The multiple-choice response item model (both unidimensional and
##' multidimensional models have the same parameterization).
##'
##' @export
##' @name Class rpf.mdim.mcm
##' @rdname rpf.mdim.mcm-class
##' @aliases rpf.mdim.mcm-class
setClass("rpf.mdim.mcm", contains='rpf.mdim',
         representation(a.prior.sdlog="numeric",
                        c.prior.alpha="numeric",
                        c.prior.beta="numeric"))

##' Convert an IRT item model name to an ID
##'
##' drm1 is the standard 3PL. drm is the multidimensional version of
##' the 3PL. gpcm1 is the Generalized Partial Credit Model.
##'
##' @param name name of the item model (string)
##' @return the integer ID assigned to the given model
rpf.id_of <- function(name) {
   .Call(get_model_names, name)
}

# similar to plink, irt.prob-class
# drm: for dichotomous response models (1PL, 2PL, and 3PL)
# gpcm: for the partial credit/generalized partial credit model
# grm: for the graded response model
# mcm: for the multiple-choice model
# nrm: for the nominal response model

setClass("rpf",
         representation(numOutcomes="numeric"))

setGeneric("rpf.prob", function(m, param, theta) standardGeneric("rpf.prob"))
setGeneric("rpf.logLik", function(m, param) standardGeneric("rpf.logLik"))
setGeneric("rpf.paramDim", function(m) standardGeneric("rpf.paramDim"))
setGeneric("rpf.rparam", function(m) standardGeneric("rpf.rparam"))
setGeneric("rpf.startingParam", function(m) standardGeneric("rpf.startingParam"))
setGeneric("rpf.getLocation", function(m,param) standardGeneric("rpf.getLocation"))
setGeneric("rpf.setLocation", function(m,param,loc) standardGeneric("rpf.setLocation"))

setClass("rpf.logistic", contains='rpf',
         representation(D="numeric"))

setClass("rpf.gpcm", contains='rpf.logistic',
         representation(a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric"))

setClass("rpf.drm", contains='rpf.logistic',
         representation(guessing="numeric",
                        a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric",
                        c.prior.alpha="numeric",
                        c.prior.beta="numeric"))

rpf.sample <- function(theta, items, params=NULL) {
    if (length(theta) == 1) {
        theta <- rnorm(theta)
    }
    numPeople <- length(theta)
    numItems <- length(items)
    if (is.null(params)) {
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

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
setGeneric("rpf.rparam", function(m) standardGeneric("rpf.rparam"))

setClass("rpf.logistic", contains='rpf',
         representation(D="numeric"))

setClass("rpf.gpcm", contains='rpf.logistic',
         representation(a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric"))

setClass("rpf.drm", contains='rpf.logistic',
         representation(a.prior.meanlog="numeric",
                        a.prior.sdlog="numeric",
                        c.prior.alpha="numeric",
                        c.prior.beta="numeric"))

rpf.sample <- function(theta, items, param=NULL) {
    if (length(theta) == 1) {
        theta <- rnorm(theta)
    }
    numPeople <- length(theta)
    numItems <- length(items)
    if (is.null(param)) {
        param <- lapply(items, rpf.rparam)
    }
    outcomes <- vapply(items, function(i) i@numOutcomes, 0)
    
    P <- array(data=0, dim=c(numItems, numPeople, max(outcomes)))
    for (ix in 1:numItems) {
        i <- items[[ix]]
        P[ix,,1:i@numOutcomes] <- rpf.prob(i, param[[ix]], theta)
    }
    ret <- apply(P, c(2,1), sample, x=1:max(outcomes), size=1, replace=F)
    return(ret)
}

# similar to plink, irt.prob-class
# drm: for dichotomous response models (1PL, 2PL, and 3PL)
# gpcm: for the partial credit/generalized partial credit model
# grm: for the graded response model
# mcm: for the multiple-choice model
# nrm: for the nominal response model

setClass("rpf",
         representation(numOutcomes="numeric"))

setGeneric("rpf.prob", function(m, param, theta) {
  standardGeneric("rpf.prob")
})

setGeneric("rpf.logLik", function(m, param) {
  standardGeneric("rpf.logLik")
})

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

#rpf.sample <- function() {}


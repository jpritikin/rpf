##' Create a dichotomous response model and associated hyperparameters.
##'
##' For discussion on the choice of priors see Cai, Yang, and
##' Hansen (2011, p. 246).
##'
##' @param numChoices the number of choices in the question
##' @param factors the number of factors
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{TRUE} when \code{factors>1} and
##' \code{FALSE} when \code{factors==1}.
##' @param a.prior.sdlog under construction
##' @param poor if TRUE, use the traditional parameterization of
##' the 1d model instead of the slope-intercept parameterization
##' @return an item model
##' @export
##' @references Cai, L., Yang, J. S., & Hansen, M. (2011). Generalized
##' Full-Information Item Bifactor Analysis.  \emph{Psychological
##' Methods, 16}(3), 221-248.
rpf.drm <- function(numChoices=5, factors=1, multidimensional, a.prior.sdlog=.5, poor=FALSE) {
  if (missing(multidimensional)) {
    multidimensional <- factors > 1
  }
  if (!multidimensional && factors > 1) {
    stop("More than 1 dimension must use a multidimensional model")
  }
  guessing <- (1/numChoices)
  m <- NULL
  id <- -1
  c.prior.sd <- .5
  c.prior.logit <- log(guessing/(1-guessing))
  if (!multidimensional) {
    id <- rpf.id_of(ifelse(poor, "drm1-", "drm1"))
    m <- new("rpf.1dim.drm",
             outcomes=2,
             factors=1,
             c.prior.logit=c.prior.logit)
  } else {
    id <- rpf.id_of("drm")
    m <- new("rpf.mdim.drm",
             outcomes=2,
             factors=factors,
             c.prior.logit=c.prior.logit)
  }
  m@spec <- c(id, 2, m@factors, a.prior.sdlog, c.prior.logit, c.prior.sd)
  m
}

### 1dim

setMethod("rpf.rparam", signature(m="rpf.1dim.drm"),
          function(m) {
            n <- 1
            c(a=rlnorm(n, meanlog=0, sdlog=.5),
              b=rnorm(n),
              c=1/(1+exp(-rnorm(n, mean=m@c.prior.logit, sd=.5))),
              u=1)
          })

### mdim

setMethod("rpf.rparam", signature(m="rpf.mdim.drm"),
          function(m) {
            c(a=rlnorm(m@factors, meanlog=0, sdlog=.5),
              b=rnorm(1),
              c=1/(1+exp(-rnorm(1, mean=m@c.prior.logit, sd=.5))),
              u=1)
          })

# Not sure if this is correct because of rotation
## as.loadings <- function(m, param) {
##   loading <- vector(mode="numeric", m@factors)
##   for (d in 1:m@factors) {
##     loading[d] <- param[d] / sqrt(1+sum(param[d:m@factors]^2))
##   }
##   loading
## }

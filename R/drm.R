##' Create a dichotomous response model and associated hyperparameters.
##'
##' For discussion on the choice of priors see Cai, Yang, and
##' Hansen (2011, p. 246).
##'
##' @param numChoices the number of choices in the question
##' @param dimensions the number of dimensions
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{TRUE} when \code{dimensions>1} and
##' \code{FALSE} when \code{dimensions==1}.
##' @return an item model
##' @export
##' @references Cai, L., Yang, J. S., & Hansen, M. (2011). Generalized
##' Full-Information Item Bifactor Analysis.  \emph{Psychological
##' Methods, 16}(3), 221-248.
rpf.drm <- function(numChoices=5, dimensions=1, multidimensional) {
  if (missing(multidimensional)) {
    multidimensional <- dimensions > 1
  }
  if (!multidimensional && dimensions > 1) {
    stop("More than 1 dimension must use a multidimensional model")
  }
  guessing <- (1/numChoices)
  m <- NULL
  id <- -1
  if (!multidimensional) {
    id <- rpf.id_of("drm1")
    m <- new("rpf.1dim.drm",
        numOutcomes=2,
        numParam=3,
        dimensions=1,
        a.prior.sdlog=.5,
        c.prior.logit=log(guessing/(1-guessing)),
        c.prior.sd=.5)
  } else {
    id <- rpf.id_of("drm")
    m <- new("rpf.mdim.drm",
        numOutcomes=2,
        dimensions=dimensions,
        numParam=2+dimensions,
        a.prior.sdlog=.5,
        c.prior.logit=log(guessing/(1-guessing)),
        c.prior.sd=.5)
  }
  m@spec <- c(id, m@numOutcomes, m@dimensions, m@a.prior.sdlog, m@c.prior.logit, m@c.prior.sd)
  m
}

### 1dim

##' @references
##' Embretson & Reise (2000, p. 184)
setMethod("rpf.info", signature(m="rpf.1dim.drm", param="numeric",
                                theta="numeric"),
          function(m, param, theta) {
            p <- rpf.prob(m, param, theta)
            a <- param[1]
            c <- param[3]
            a^2 * (p[,1]/p[,2]) * ((p[,2]-c)^2/(1-c)^2)
          })

setMethod("rpf.rparam", signature(m="rpf.1dim.drm"),
          function(m) {
            n <- 1
            c(a=rlnorm(n, meanlog=0, sdlog=m@a.prior.sdlog),
              b=rnorm(n),
              c=1/(1+exp(-rnorm(n, mean=m@c.prior.logit, sd=m@c.prior.sd))))
          })

### mdim

setMethod("rpf.rparam", signature(m="rpf.mdim.drm"),
          function(m) {
            c(a=rlnorm(m@dimensions, meanlog=0, sdlog=m@a.prior.sdlog),
              b=rnorm(1),
              c=1/(1+exp(-rnorm(1, mean=m@c.prior.logit, sd=m@c.prior.sd))))
          })

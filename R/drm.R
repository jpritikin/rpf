##' Create a dichotomous response model
##'
##' For discussion on the choice of priors see Cai, Yang, and
##' Hansen (2011, p. 246).
##'
##' @param factors the number of factors
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{TRUE} when \code{factors>1} and
##' \code{FALSE} when \code{factors==1}.
##' @param poor if TRUE, use the traditional parameterization of
##' the 1d model instead of the slope-intercept parameterization
##' @return an item model
##' @export
##' @references Cai, L., Yang, J. S., & Hansen, M. (2011). Generalized
##' Full-Information Item Bifactor Analysis.  \emph{Psychological
##' Methods, 16}(3), 221-248.
rpf.drm <- function(factors=1, multidimensional=TRUE, poor=FALSE) {
  if (!multidimensional && factors > 1) {
    stop("More than 1 dimension must use a multidimensional model")
  }
  m <- NULL
  id <- -1
  if (!multidimensional) {
    if (!poor) stop("The old parameterization is no longer available")
    id <- rpf.id_of("drm1-")
    m <- new("rpf.1dim.drm",
             outcomes=2,
             factors=1)
  } else {
    id <- rpf.id_of("drm")
    m <- new("rpf.mdim.drm",
             outcomes=2,
             factors=factors)
  }
  m@spec <- c(id, 2, m@factors)
  m
}

### 1dim

setMethod("rpf.rparam", signature(m="rpf.1dim.drm"),
          function(m) {
            n <- 1
            c(a=rlnorm(n, meanlog=0, sdlog=.5),
              b=rnorm(n),
              c=rbeta(1, 5,17),
              u=rbeta(1, 17,5))
          })

### mdim

setMethod("rpf.rparam", signature(m="rpf.mdim.drm"),
          function(m) {
            c(a=rlnorm(m@factors, meanlog=0, sdlog=.5),
              b=rnorm(1),
              c=rbeta(1, 5,17),
              u=rbeta(1, 17,5))
          })

# Not sure if this is correct because of rotation
## as.loadings <- function(m, param) {
##   loading <- vector(mode="numeric", m@factors)
##   for (d in 1:m@factors) {
##     loading[d] <- param[d] / sqrt(1+sum(param[d:m@factors]^2))
##   }
##   loading
## }

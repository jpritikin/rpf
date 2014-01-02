##' Create a graded response model
##'
##' For outcomes k in 0 to K, slope vector a, intercept vector c, and latent ability vector theta,
##' the response probability function is
##' \deqn{\mathrm P(\mathrm{pick}=0|a,c,\theta) = 1- \mathrm P(\mathrm{pick}=1|a,c_1,\theta)
##' }{P(pick=0|a,c,th) = 1-P(pick=1|a,c_1,th)}
##' \deqn{\mathrm P(\mathrm{pick}=k|a,c,\theta) = \frac{1}{1+\exp(-(a\theta + c_k))} - \frac{1}{1+\exp(-(a\theta + c_{k+1}))}
##' }{P(pick=k|a,c,th) = 1/(1+exp(-(a th + c_k))) - 1/(1+exp(-(a th + c_(k+1))))}
##' \deqn{\mathrm P(\mathrm{pick}=K|a,c,\theta) = \frac{1}{1+\exp(-(a\theta + c_K))}
##' }{P(pick=K|a,c,th) = 1/(1+exp(-(a th + c_K)))}
##'
##' The graded response model was designed for a item with a series of
##' dependent parts where a higher score implies that easier parts of
##' the item were surmounted. If there is any chance your polytomous
##' item has independent parts then consider \code{\link{rpf.nrm}}.
##' If your categories cannot cross then the graded response model
##' provides a little more information than the nominal model.
##' Stronger a priori assumptions offer provide more power at the cost
##' of flexibility.
##' 
##' @param outcomes The number of choices available
##' @param factors the number of factors
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{TRUE}.
##' @return an item model
##' @export
rpf.grm <- function(outcomes=2, factors=1, multidimensional=TRUE) {
  if (!multidimensional && factors > 1) {
    stop("More than 1 dimension must use a multidimensional model")
  }
  m <- NULL
  id <- -1
  if (!multidimensional) {
    stop("The old parameterization is no longer available")
  } else {
    id <- rpf.id_of("grm")
    m <- new("rpf.mdim.grm",
             outcomes=outcomes,
             factors=factors)
  }
  m@spec <- c(id, m@outcomes, m@factors)
  m
}

### mdim

setMethod("rpf.modify", signature(m="rpf.mdim.graded", factors="numeric"),
          function(m, factors) {
              rpf.grm(m@outcomes, factors)
          })

setMethod("rpf.rparam", signature(m="rpf.mdim.graded"),
          function(m) {
              a <- rlnorm(m@factors, meanlog=0, sdlog=.5)
              b <- rnorm(m@outcomes-1)
              b <- b[order(-b)]
              c(a=a,b=b)
          })

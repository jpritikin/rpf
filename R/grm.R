##' Create a graded response model
##'
##' This function instantiates a graded response model.
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
##' Defaults to \code{TRUE} when \code{factors>1} and
##' \code{FALSE} when \code{factors==1}.
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

### 1dim

setMethod("rpf.rparam", signature(m="rpf.1dim.graded"),
          function(m) {
              a <- rlnorm(1, meanlog=0, sdlog=.5)
              b <- sort(rnorm(m@outcomes-1))
              c(a=a,b=b)
          })

### mdim

setMethod("rpf.rparam", signature(m="rpf.mdim.graded"),
          function(m) {
              a <- rlnorm(m@factors, meanlog=0, sdlog=.5)
              b <- rnorm(m@outcomes-1)
              b <- b[order(-b)]
              c(a=a,b=b)
          })

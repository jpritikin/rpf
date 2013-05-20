##' Create a graded response model and associated hyperparameters.
##'
##' This function instantiates a graded response model. Bayesian
##' priors are only used to generate plausible random parameters.
##'
##' The graded response model was designed for a item with a series of
##' dependent parts where a higher score implies that easier parts of
##' the item were surmounted. If there is any chance your polytomous
##' item has independent parts then consider \code{\link{rpf.gpcm}}.
##' If your categories cannot cross then the graded response model
##' provides a little more information than the GPCM. Stronger a
##' priori assumptions offer provide more power at the cost of
##' flexibility.
##' 
##' @param numOutcomes The number of choices available
##' @param dimensions the number of dimensions
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{TRUE} when \code{dimensions>1} and
##' \code{FALSE} when \code{dimensions==1}.
##' @return an item model
##' @export
rpf.grm <- function(numOutcomes=2, dimensions=1, multidimensional) {
  if (missing(multidimensional)) {
    multidimensional <- dimensions > 1
  }
  if (!multidimensional && dimensions > 1) {
    stop("More than 1 dimension must use a multidimensional model")
  }
  m <- NULL
  id <- -1
  if (!multidimensional) {
    id <- rpf.id_of("grm1")
    m <- new("rpf.1dim.grm",
             numOutcomes=numOutcomes,
             dimensions=1,
             a.prior.sdlog=.5)
  } else {
    id <- rpf.id_of("grm")
    m <- new("rpf.mdim.grm",
             numOutcomes=numOutcomes,
             dimensions=dimensions,
             a.prior.sdlog=.5)
  }
  m@spec <- c(id, m@numOutcomes, m@dimensions)
  m
}


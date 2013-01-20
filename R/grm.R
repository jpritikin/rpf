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
##' @author Jonathan Weeks <weeksjp@@gmail.com>
rpf.grm <- function(numOutcomes=2, dimensions=1, multidimensional) {
  if (missing(multidimensional)) {
    multidimensional <- dimensions > 1
  }
  if (!multidimensional && dimensions > 1) {
    stop("More than 1 dimension must use a multidimensional model")
  }
  if (!multidimensional) {
    new("rpf.1dim.grm",
        numOutcomes=numOutcomes,
        dimensions=1,
        numParam=numOutcomes,
        a.prior.sdlog=.5)
  } else {
    new("rpf.mdim.grm",
        numOutcomes=numOutcomes,
        dimensions=dimensions,
        numParam=dimensions + numOutcomes - 1,
        a.prior.sdlog=.5)
  }
}

### 1dim

setMethod("rpf.prob", signature(m="rpf.1dim.grm", param="numeric",
                                theta="numeric"),
          function(m, param, theta) {
            a <- param[1]
            b <- param[-1]

            ct <- m@numOutcomes - 1

            p <- array(dim=c(length(theta), m@numOutcomes))
            
            p[,1] <- 1-1/(1+exp(-a*(theta-b[1])))
            
            for (k in 1:ct) {
              if (k<ct) {
                cp <- (1/(1+exp(-a*(theta-b[k])))) - (1/(1+exp(-a*(theta-b[k+1]))))
              } else if (k==ct) {
                cp <- 1/(1+exp(-a*(theta-b[k])))
              }
              p[,k+1] <- cp
            }
            return(p)
          })

### mdim

setMethod("rpf.prob", signature(m="rpf.mdim.grm", param="numeric",
                                theta="matrix"),
          function(m, param, theta) {
            a <- param[1:m@dimensions]
            b <- param[-1:-m@dimensions]

            ct <- m@numOutcomes - 1
            
            numPersons <- dim(theta)[1]
            p <- array(dim=c(numPersons, m@numOutcomes))

            ##   Compute the probabilities for the lowest category
            p[,1] <- 1-1/(1+exp(-(theta %*% a+b[1])))
            
            for (k in 1:ct) {
              if (k<ct) {
                cp <- (1/(1+exp(-(theta %*% a+b[k]))))-(1/(1+exp(-(theta %*% a+b[k+1]))))
              } else if (k==ct) {
                ##   Compute the probabilities for the highest category
                cp <- 1/(1+exp(-(theta %*% a+b[k])))
              }
              p[,k+1] <- cp
            }
            return(p)
          })

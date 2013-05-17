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
            nc <- m@numOutcomes
            a <- param[1:m@dimensions]
            c <- param[-m@dimensions:-1]
            
            numPersons <- dim(theta)[1]
            
            T <- matrix(0,numPersons,nc+1)
            T[,1] <- 1
            T[,nc+1] <- 0
            for (k in 1:(nc-1)) {
              z <- c[k] + a %*% t(theta)
              T[,k+1] <- 1/(1+exp(-z))
            }
            P <- matrix(0,numPersons,nc)
            for (k in 1:nc) {
              P[,k] <- T[,k]-T[,k+1]
            }
            P
          })

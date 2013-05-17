# Copied from flexmirt.R
Tnom.trend <- function(nc) {
  T <- matrix(0,nc,nc-1)
  for (i in 1:nc) {
    T[i,1] <- i-1
  }
  for (k in 2:(nc-1)) {
    for (i in 2:(nc-1)) {
      T[k,i] <- sin(pi*(i-1)*(k-1)/(nc-1))
    }
  }
  return(T)
}

# Copied from flexmirt.R
Tnom.id <- function(nc) {
  T <- matrix(0,nc,nc-1)
  T[nc,1] <- nc-1
  T[2:(nc-1),2:(nc-1)] <- diag(nc-2)
  return(T)
}

build.T <- function(numOutcomes, got) {
  if (!is.matrix(got)) {
    if (got == "id") {
      got <- Tnom.id(numOutcomes)
    } else if (got == "trend") {
      got <- Tnom.trend(numOutcomes)
    } else {
      stop(paste("T matrix", deparse(got), "not recognized"))
    }
  }
  if (all(dim(got) == c(numOutcomes,numOutcomes-1))) {
    if (any(got[1,] != 0)) warn("Non-zero T[1,] will be ignored")
    got <- got[-1,]
  }
  if (all(dim(got) != rep(numOutcomes-1, 2))) {
    stop(paste("T matrix must be of dimensions",
               paste(rep(numOutcomes-1, 2), collapse="x"),
               "not", paste(dim(got), collapse="x")))
  }
  got
}

##' Create a nominal response model and associated hyperparameters.
##'
##' This function instantiates a nominal response model. Bayesian
##' priors are only used to generate plausible random parameters.
##' 
##' @param numOutcomes The number of choices available
##' @param dimensions the number of dimensions
##' @return an item model
##' @export
rpf.nrm <- function(numOutcomes=3, dimensions=1, T.a="trend", T.c="trend") {
  T.a <- build.T(numOutcomes, T.a)
  T.c <- build.T(numOutcomes, T.c)
  id <- rpf.id_of("nominal")
  m <- new("rpf.mdim.nrm",
           numOutcomes=numOutcomes,
           dimensions=dimensions,
           a.prior.sdlog=.5)
  m@spec <- c(id, numOutcomes, dimensions, T.a, T.c)
  m
}

getT <- function(m, tx) {
  Tsize <- (m@numOutcomes-1L)^2
  offset <- 4 + tx * Tsize
  matrix(m@spec[offset:(offset+Tsize-1)], m@numOutcomes-1L, m@numOutcomes-1L)
}

setMethod("rpf.rparam", signature(m="rpf.mdim.nrm"),
          function(m) {
            a <- rlnorm(m@dimensions, sdlog=.5)
            ak <- abs(rnorm(m@numOutcomes-1, mean=1, sd=.25))
            ck <- sort(rnorm(m@numOutcomes-1))
            c(a=a,
              alf=solve(getT(m,0)) %*% ak,
              gam=solve(getT(m,1)) %*% ck)
          })

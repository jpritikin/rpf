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
  return(T[-1,])
}

# Copied from flexmirt.R
Tnom.id <- function(nc) {
  T <- matrix(0,nc,nc-1)
  T[nc,1] <- nc-1
  T[2:(nc-1),2:(nc-1)] <- diag(nc-2)
  return(T[-1,])
}

build.T <- function(outcomes, got) {
  if (!is.matrix(got)) {
    if (got == "id") {
      got <- Tnom.id(outcomes)
    } else if (got == "trend") {
      got <- Tnom.trend(outcomes)
    } else if (got == "random") {
      while (1) {
        side <- outcomes-1
        got <- matrix(rnorm(side*side), side, side)
        invertible <- try(solve(got), silent=TRUE)
        if (!inherits(invertible, "try-error")) break
      }
    } else {
      stop(paste("T matrix", deparse(got), "not recognized"))
    }
  }
  if (all(dim(got) == c(outcomes,outcomes-1))) {
    if (any(got[1,] != 0)) warning("Non-zero T[1,] will be ignored")
    got <- got[-1,]
  }
  if (all(dim(got) != rep(outcomes-1, 2))) {
    stop(paste("T matrix must be of dimension",
               paste(rep(outcomes-1, 2), collapse="x"),
               "not", paste(dim(got), collapse="x")))
  }
  got
}

##' Create a nominal response model
##'
##' This function instantiates a nominal response model. The T matrix
##' must be an invertible square matrix of dimension outcomes-1. As a
##' shortcut, either T matrix can be specified as "trend" for a
##' Fourier basis or as "id" for an identity basis.
##' 
##' @param outcomes The number of choices available
##' @param factors the number of factors
##' @param T.a the T matrix for slope parameters
##' @param T.c the T matrix for intercept parameters
##' @return an item model
##' @references Thissen, D., Cai, L., & Bock, R. D. (2010). The
##' Nominal Categories Item Response Model. In M. L. Nering &
##' R. Ostini (Eds.), \emph{Handbook of Polytomous Item Response
##' Theory Models} (pp. 43--75). Routledge.
rpf.nrm <- function(outcomes=3, factors=1, T.a="trend", T.c="trend") {
  T.a <- build.T(outcomes, T.a)
  T.c <- build.T(outcomes, T.c)
  id <- rpf.id_of("nominal")
  m <- new("rpf.mdim.nrm",
           outcomes=outcomes,
           factors=factors)
  m@spec <- c(id, outcomes, factors, T.a, T.c, solve(T.a), solve(T.c))
  m
}

getT <- function(m, tx) {
  Tsize <- (m@outcomes-1L)^2
  offset <- 4 + tx * Tsize
  matrix(m@spec[offset:(offset+Tsize-1)], m@outcomes-1L, m@outcomes-1L)
}

setMethod("rpf.rparam", signature(m="rpf.mdim.nrm"),
          function(m) {
            a <- rlnorm(m@factors, sdlog=.5)
            ak <- abs(rnorm(m@outcomes-1, mean=1, sd=.25))
            ck <- sort(rnorm(m@outcomes-1))
            c(a=a,
              alf=getT(m,2) %*% ak,
              gam=getT(m,3) %*% ck)
          })

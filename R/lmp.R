##' Create logistic function of a monotonic polynomial (LMP) model
##'
##' This model is a dichotomous response model originally proposed by
##' Liang (2007) and is implemented using the parameterization by
##' Falk & Cai (in press).
##'
##' The LMP model replaces the linear predictor part of the
##' two-parameter logistic function with a monotonic polynomial,
##' \eqn{m(\theta,\omega,\xi,\mathbf{\alpha},\mathbf{\tau})},
##'
##' \deqn{\mathrm P(\mathrm{pick}=1|\omega,\xi,\mathbf{\alpha},\mathbf{\tau},\theta)
##' = \frac{1}{1+\exp(-(\xi + m(\theta,\omega,\xi,\mathbf{\alpha},\mathbf{\tau})))}
##' }{P(pick=1|omega,xi,alpha,tau,th) = 1/(1+exp(-(xi + m(theta,omega,xi,alpha,tau))))}
##'
##' where \eqn{\mathbf{\alpha}}{alpha} and \eqn{\mathbf{\tau}}{tau} are vectors
##' of length k. See Falk & Cai (in press) for more details as to how the
##' polynomial is constructed.
##'
##' The order of the polynomial is always odd and is controlled by
##' the user specified non-negative integer, k. The model contains
##' 2+2*k parameters and are used as input to the \code{rpf.prob}
##' function in the following order:
##' \eqn{\omega}{omega} - the natural log of the slope of the item model when k=0,
##' \eqn{\xi}{xi} - the intercept,
##' \eqn{\alpha}{alpha} and \eqn{\tau}{tau} - two parameters that control bends in
##' the polynomial. These latter parameters are repeated in the same order for
##' models with k>1.
##'
##' At the lowest order polynomial (k=0) the model reduces to the
##' two-parameter logistic model. Due to this, k=0 may not currently be
##' supported.
##'
##' @param k a non-negative integer that controls the order of the
##' polynomial (2k*1) with a default of k=1 (3rd order polynomial).
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{FALSE}. The multidimensional version is not yet
##' available.
##' @return an item model
##' @export
##' @references Falk, C. F., & Cai, L. (in press). Maximum marginal likelihood
##' estimation of a monotonic polynomial generalized partial credit model with
##' applications to multiple group analysis. \emph{Psychometrika}.
##'
##' Liang (2007). \emph{A semi-parametric approach to estimating item response
##' functions}. Unpublished doctoral dissertation, Department of Psychology,
##' The Ohio State University.
##' @examples
##' spec <- rpf.lmp(1) # 3rd order polynomial
##' theta<-seq(-3,3,.1)
##' p<-rpf.prob(spec, c(-.11,.37,.24,-.21),theta)
##'
##' spec <- rpf.lmp(2) # 5th order polynomial
##' p<-rpf.prob(spec, c(.69,.71,-.5,-8.48,.52,-3.32),theta)

## My own version
rpf.lmp <- function(k=1, multidimensional=FALSE) {
  if(!(k%%1==0)|!(k>=1)){
    stop("k must be an integer > 0")
  }
  if(multidimensional){
      stop("Multidimensional LMPA model is not yet supported")
  }
  m <- NULL
  id <- -1
  id <- rpf.id_of("lmp")
  m <- new("rpf.lmp.drm",
           outcomes=2,
           factors=1)
  m@spec <- c(id, 2, m@factors, k)
  m
}

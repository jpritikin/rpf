##' Create monotonic polynomial graded response (GR-MP) model
##'
##' The GR-MP model replaces the linear predictor of the graded
##' response model with a monotonic polynomial. In contrast to the
##' models by Falk and Cai (2016), the GR-MP was implmented using
##' a slightly different parameterization for the polynomial. In 
##' particular, the polynomial is parameterized such that boundary
##' descrimination functions for the GR-MP will be all monotonically
##' increasing or decreasing for any given item. This allows the
##' possibility of items that load either negatively or positively
##' on the latent trait.
##'
##' The order of the polynomial is always odd and is controlled by
##' the user specified non-negative integer, k. The model contains
##' 1+(outcomtes-1)+2*k parameters and are used as input to the \code{\link{rpf.prob}}
##' function in the following order:
##' \eqn{\lambda}{lambda} - slope of the item model when k=0,
##' \eqn{\xi}{xi} - a (outcomes-1)-length vector of intercept parameters,
##' \eqn{\alpha}{alpha} and \eqn{\tau}{tau} - two parameters that control bends in
##' the polynomial. These latter parameters are repeated in the same order for
##' models with k>1.  For example, a k=2 polynomial with 3 categories will have an item
##' parameter vector of: \eqn{\lambda, \xi_1, \xi_2, \alpha_1, \tau_1, \alpha_2, \tau_2}{
##' lambda, xi1, xi2, alpha1, tau1, alpha2, tau2}.
##'
##'
##' @param outcomes The number of possible response categories.
##' @param k a non-negative integer that controls the order of the
##' polynomial (2k+1) with a default of k=0 (1st order polynomial = generalized partial credit model).
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{FALSE}. The multidimensional version is not yet
##' available.
##' @return an item model
##' @references Falk, C. F., & Cai, L. (2016). Maximum marginal likelihood
##' estimation of a monotonic polynomial generalized partial credit model with
##' applications to multiple group analysis. \emph{Psychometrika, 81}, 434-460.
##' \url{http://dx.doi.org/10.1007/s11336-014-9428-7}
##' @seealso \code{\link{rpf.lmp}} \code{\link{rpf.gpcmp}}
##'
##' @examples
##' spec <- rpf.grmp(5,2) # 5-category, 3rd order polynomial
##' theta<-seq(-3,3,.1)
##' p<-rpf.prob(spec, c(2.77,2,1,0,-1,.89,-8.7,-.74,-8.99),theta)

rpf.grmp <- function(outcomes=2, k=0, multidimensional=FALSE) {
  if(!(k%%1==0)){
    stop("k must be an integer >= 0")
  }
  if(multidimensional){
      stop("Multidimensional grmp model is not yet supported")
  }
  m <- NULL
  id <- -1
  id <- rpf.id_of("grmp")
  m <- new("rpf.1dim.grmp",
           outcomes=outcomes,
           factors=1)
  m@spec <- c(id, m@outcomes, m@factors, k)
  m
}

setMethod("rpf.rparam", signature(m="rpf.1dim.grmp"),
          function(m, version) {
            n <- 1
            k<-m$spec[4] ## ok to hardcode this index?
            ret<-c(lambda=rlnorm(n, 0, .5)) # random overall slope
            # randomly generate xi
            xi <- rnorm(m@outcomes-1)
            xi <- xi[order(-xi)]
            ret<-c(ret,xi=xi)
            if(k>0){
                for(i in 1:k){
                    ret<-c(ret,runif(n,-1,1),log(runif(n,.0001,1)))
                    names(ret)[(m@outcomes+1+(i-1)*2):(m@outcomes+(i*2))]<-c(paste("alpha",i,sep=""),paste("tau",i,sep=""))
                }
            }
            ret
        })

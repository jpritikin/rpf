##' Compute the sum-score EAP table
##'
##' TODO: Optimize for two-tier covariance structure
##' 
##' @param grp a list with spec, param, mean, and cov
##' @param width positive width of quadrature in Z units
##' @param pts number of quadrature points
##' @examples
##' # see Thissen, Pommerich, Billeaud, & Williams (1995, Table 2)
##'  spec <- list()
##'  spec[1:3] <- rpf.grm(outcomes=4)
##'  
##'  param <- matrix(c(1.87, .65, 1.97, 3.14,
##'                    2.66, .12, 1.57, 2.69,
##'                    1.24, .08, 2.03, 4.3), nrow=4)
##'  # fix parameterization
##'  param <- apply(param, 2, function(p) c(p[1], p[2:4] * -p[1]))
##'  
##'  grp <- list(spec=spec, mean=0, cov=matrix(1,1,1), param=param)
##'  sumScoreEAP(grp)
sumScoreEAP <- function(grp, width=6.0, pts=49L) {
	.Call(ssEAP_wrapper, grp, width, pts)
}

##' Compute the observed sum-score
##'
##' @param grp a list with spec, param, and data
##' @param mask a vector of logicals indicating which items to include
##' @examples
##' spec <- list()
##' spec[1:3] <- rpf.grm(outcomes=3)
##' param <- sapply(spec, rpf.rparam)
##' data <- rpf.sample(5, spec, param)
##' colnames(param) <- colnames(data)
##' grp <- list(spec=spec, param=param, data=data)
##' observedSumScore(grp, rep(TRUE, length(spec)))
observedSumScore <- function(grp, mask) {
	.Call(observedSumScore_wrapper, grp, mask)
}

##' Produce an item outcome by observed sum-score table
##'
##' @param grp a list with spec, param, and data
##' @param mask a vector of logicals indicating which items to include
##' @param interest index or name of the item of interest
##' @examples
##' set.seed(1)
##' spec <- list()
##' spec[1:3] <- rpf.grm(outcomes=3)
##' param <- sapply(spec, rpf.rparam)
##' data <- rpf.sample(5, spec, param)
##' colnames(param) <- colnames(data)
##' grp <- list(spec=spec, param=param, data=data)
##' itemOutcomeBySumScore(grp, c(FALSE,TRUE,TRUE), 1L)
itemOutcomeBySumScore <- function(grp, mask, interest) {
	if (is.character(interest)) {
		interest <- match(interest, colnames(grp$param))
	}
	tbl <- .Call(itemOutcomeBySumScore_wrapper, grp, mask, interest)
	rownames(tbl) <- 0:(nrow(tbl)-1L)
	tbl
}

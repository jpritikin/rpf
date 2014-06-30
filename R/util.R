ssEAP <- function(grp, qwidth, qpoints, mask, twotier=FALSE, debug=FALSE) {
	if (missing(mask)) {
		mask <- rep(TRUE, ncol(grp$param))
	}
	.Call(ssEAP_wrapper, grp, qwidth, qpoints, mask, twotier, debug)
}

##' Compute the sum-score EAP table
##'
##' Observed tables cannot be computed when data is
##' missing. Therefore, you can optionally omit items with the
##' greatest number of responses missing when conducting the
##' distribution test.
##' 
##' When two-tier covariance structure is detected, EAP scores are
##' only reported for primary factors. It is possible to compute EAP
##' scores for specific factors, but it is not clear why this would be
##' useful because they are conditional on the specific factor sum
##' scores. Moveover, the algorithm to compute them efficiently has not been
##' published yet (as of Jun 2014).
##'
##' @param grp a list with spec, param, mean, and cov
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param qwidth positive width of quadrature in Z units
##' @param qpoints number of quadrature points
##' @param distributionTest whether to perform the latent distribution test
##' @param omit number of items to omit from the latent distribution test
##' @param .twotier whether to enable the two-tier optimization
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
##' @references
##' Li, Z., & Cai, L. (2012, July). Summed score likelihood based indices for testing
##' latent variable distribution fit in Item Response Theory. Paper presented at
##' the annual International Meeting of the Psychometric Society, Lincoln,
##' NE. Retrieved from http://www.cse.ucla.edu/downloads/files/SD2-final-4.pdf
sumScoreEAP <- function(grp, ..., qwidth=6.0, qpoints=49L, distributionTest=NULL, omit=0L, .twotier=TRUE) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}

	if (missing(qwidth) && !is.null(grp$qwidth)) { qwidth <- grp$qwidth }
	if (missing(qpoints) && !is.null(grp$qpoints)) { qpoints <- grp$qpoints }

	tbl <- ssEAP(grp, qwidth, qpoints, twotier=.twotier)
	rownames(tbl) <- 0:(nrow(tbl)-1)
	result <- list(tbl=tbl, distributionTest=FALSE)
	if ((is.null(distributionTest) && !is.null(grp$data)) || (!is.null(distributionTest) && distributionTest)) {
		if (!is.null(distributionTest) && is.null(grp$data)) {
			stop("distributionTest cannot be conducted because there is no data")
		}
		result$distributionTest <- TRUE
		mask <- rep(TRUE, ncol(grp$param))
		if (omit > 0) {
			if (omit >= ncol(grp$param)) stop("Cannot omit all the items")
			nacount <- sort(-sapply(grp$data, function(x) sum(is.na(x))))
			omit <- min(sum(nacount != 0), omit)
			toOmit <- names(nacount)[1:omit]
			mask[match(toOmit, colnames(grp$param))] <- FALSE
			tbl <- ssEAP(grp, qwidth, qpoints, mask)
			result$omitted <- toOmit
		}
		oss <- observedSumScore(grp, mask=mask)
		result$n <- oss$n
		result$observed <- oss$dist
		obs <- matrix(oss$dist, ncol=1)
		size <- sum(obs)
		expected <- matrix(size * tbl[,1], ncol=1)
		result$rms.p <- log(ptw2011.gof.test(obs, expected))

		kc <- .Call(collapse_wrapper, obs, expected)
		obs <- kc$O
		expected <- kc$E
		mask <- !is.na(expected) & expected!=0
		result$pearson.chisq <- sum((obs[mask] - expected[mask])^2 / expected[mask])
		result$pearson.df <- sum(mask)-1L
		result$pearson.p <- pchisq(result$pearson.chisq, result$pearson.df, lower.tail=FALSE, log.p=TRUE)
	}
	class(result) <- "summary.sumScoreEAP"
	result
}

print.summary.sumScoreEAP <- function(x,...) {
	print(x$tbl)
	if (x$distributionTest) {
		cat(sprintf("\nLatent distribution fit test (n=%d):\n", x$n))
	}
	if (!is.null(x$omitted)) {
		cat(paste("  Omitted:", paste(x$omitted, collapse=" "), "\n"))
	}
	if (!is.null(x$rms.p)) {
		cat(sprintf("  RMS log(p) = %.2f\n", x$rms.p))
	}
	if (!is.null(x$pearson.p)) {
		cat(sprintf("  Pearson X^2(%3d) = %.2f, log(p) = %.2f\n",
			    x$pearson.df, x$pearson.chisq, x$pearson.p))
	}
}

##' Compute the observed sum-score
##'
##' @param grp a list with spec, param, and data
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param mask a vector of logicals indicating which items to include
##' @param summary whether to return a summary (default) or per-row scores
##' @examples
##' spec <- list()
##' spec[1:3] <- rpf.grm(outcomes=3)
##' param <- sapply(spec, rpf.rparam)
##' data <- rpf.sample(5, spec, param)
##' colnames(param) <- colnames(data)
##' grp <- list(spec=spec, param=param, data=data)
##' observedSumScore(grp)
observedSumScore <- function(grp, ..., mask, summary=TRUE) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}
	if (missing(mask)) {
		mask <- rep(TRUE, ncol(grp$param))
	}
	if (!summary) {
		cols <- colnames(grp$param)[mask]
		dat <- grp$data[,cols]
		ss <- apply(sapply(dat, unclass) - 1, 1, sum)
		names(ss) <- rownames(dat)
		return(ss)
	}
	got <- .Call(observedSumScore_wrapper, grp, mask)
	class(got) <- "summary.observedSumScore"
	got
}

print.summary.observedSumScore <- function(x,...) {
	print(x$dist)
	cat(sprintf("  N = %d\n", x$n))
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
	got <- .Call(itemOutcomeBySumScore_wrapper, grp, mask, interest)
	rownames(got$table) <- 0:(nrow(got$table)-1L)
	col <- colnames(grp$param)[interest]
	colnames(got$table) <- levels(grp$data[,col])
	class(got) <- "summary.itemOutcomeBySumScore"
	got
}

print.summary.itemOutcomeBySumScore <- function(x,...) {
	print(x$table)
	cat(sprintf("  N = %d\n", x$n))
}

##' Compute EAP scores
##'
##' @param grp a list with spec, param, and data
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param naAction action for rows with fewer than
##' \code{minItemsPerScore}. Defaults to 'fail'. If 'pass', will fill
##' with NAs.
##' @examples
##' spec <- list()
##' spec[1:3] <- rpf.grm(outcomes=3)
##' param <- sapply(spec, rpf.rparam)
##' data <- rpf.sample(5, spec, param)
##' colnames(param) <- colnames(data)
##' grp <- list(spec=spec, param=param, data=data, minItemsPerScore=1L)
##' EAPscores(grp)
EAPscores <- function(grp, ..., naAction="fail") {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}

	validAction <- c('pass', 'fail')
	if (is.na(match(naAction, validAction))) {
		stop(paste("naAction must be one of", paste(validAction, collapse=", ")))
	}

	.Call(eap_wrapper, grp, naAction == "fail")
}

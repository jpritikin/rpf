##' Identify the columns with most missing data
##'
##' If a reference column is given then only rows that are not missing
##' on the reference column are considered. Otherwise all rows are
##' considered.
##' 
##' @param grp an IFA group
##' @param omit the maximum number of items to omit
##' @param ref the reference column (optional)
bestToOmit <- function(grp, omit, ref=NULL) {
	if (missing(omit)) stop("How many items to omit?")
	if (omit == 0) return(NULL)
	dat <- grp$data
	wcol <- 1
	if (!is.null(grp$weightColumn)) {
		wcol <- dat[[grp$weightColumn]]
		dat <- dat[,-match(grp$weightColumn, colnames(dat))]
	}
	if (omit >= ncol(dat)) stop("Cannot omit all columns")
	if (!is.null(ref)) {
		dat <- dat[!is.na(dat[[ref]]),]
	}
	nacount <- apply(dat, 2, function(c) sum(is.na(c) * wcol))
	omit <- min(omit, sum(nacount > 0))
	if (omit == 0) return(NULL)
	names(sort(-nacount)[1:omit])
}

##' Omit the given items
##'
##' @param grp an IFA group
##' @param excol vector of column names to omit
omitItems <- function(grp, excol) {
	if (missing(excol)) stop("Which items to omit?")
	if (length(excol) == 0) return(grp)
	imask <- -match(excol, colnames(grp$param))
	grp$spec <- grp$spec[imask]
	grp$param <- grp$param[,imask]
	grp$free <- grp$free[,imask]
	grp$data <- grp$data[,-match(excol, colnames(grp$data))]

	# We need to repack the data because
	# rows that only differed on the column
	# we removed are now the same.
	if (!is.null(grp$weightColumn)) {
		data <- expandDataFrame(grp$data, grp$weightColumn)
		grp$data <- compressDataFrame(data)
	}

	if (!is.null(grp$observedStats)) {
		grp$observedStats <- nrow(grp$data)
	}
	grp$omitted <- c(grp$omitted, excol)
	grp
}

##' Omit items with the most missing data
##'
##' Items with no missing data are never omitted, regardless of the
##' number of items requested.
##' 
##' @param grp an IFA group
##' @param omit the maximum number of items to omit
omitMostMissing <- function(grp, omit) {
	omitItems(grp, bestToOmit(grp, omit))
}

ssEAP <- function(grp, qwidth, qpoints, mask, twotier=FALSE, debug=FALSE) {
	if (missing(mask)) {
		mask <- rep(TRUE, ncol(grp$param))
	}
	.Call(ssEAP_wrapper, grp, qwidth, qpoints, mask, twotier, debug)
}

sumScoreEAPTestInternal <- function(result) {
	expected <- matrix(result$expected, ncol=1)
	obs <- matrix(result$observed, ncol=1)

	result$rms.p <- log(ptw2011.gof.test(obs, expected))

	kc <- .Call(collapse_wrapper, obs, expected)
	obs <- kc$O
	expected <- kc$E
	mask <- !is.na(expected) & expected!=0
	result$pearson.chisq <- sum((obs[mask] - expected[mask])^2 / expected[mask])
	result$pearson.df <- sum(mask)-1L
	result$pearson.p <- pchisq(result$pearson.chisq, result$pearson.df, lower.tail=FALSE, log.p=TRUE)

	class(result) <- "summary.sumScoreEAPTest"
	result
}

##' Conduct the sum-score EAP distribution test
##' 
##' @param grp a list with spec, param, mean, and cov
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param .twotier whether to enable the two-tier optimization
##' @references
##' Li, Z., & Cai, L. (2012, July). Summed score likelihood based indices for testing
##' latent variable distribution fit in Item Response Theory. Paper presented at
##' the annual International Meeting of the Psychometric Society, Lincoln,
##' NE. Retrieved from http://www.cse.ucla.edu/downloads/files/SD2-final-4.pdf
sumScoreEAPTest <- function(grp, ..., .twotier=TRUE) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}
	if (is.null(grp$data)) {
		stop("distributionTest cannot be conducted because there is no data")
	}
	qwidth <- grp$qwidth
	if (is.null(qwidth)) { qwidth <- 6 }
	qpoints <- grp$qpoints
	if (is.null(qpoints)) { qpoints <- 49 }
	tbl <- ssEAP(grp, qwidth, qpoints, twotier=.twotier)
	rownames(tbl) <- 0:(nrow(tbl)-1)
	result <- list(tbl=tbl)
	oss <- observedSumScore(grp)
	result$n <- oss$n
	result$observed <- oss$dist
	result$expected <- result$n * tbl[,1]
	names(result$observed) <- rownames(tbl)
	names(result$expected) <- rownames(tbl)
	result$omitted <- grp$omitted
	result <- sumScoreEAPTestInternal(result)
	result
}

"+.summary.sumScoreEAPTest" <- function(e1, e2) {
	e2name <- deparse(substitute(e2))
	if (!inherits(e2, "summary.sumScoreEAPTest")) {
		stop("Don't know how to add ", e2name, " to a sumScoreEAPTest",
		     call. = FALSE)
	}

	if (length(e1$observed) != length(e2$observed)) {
		stop("The two groups have a different maximum sum-score. Sum-score tests cannot be combined")
	}
	if (any(e1$omitted != e2$omitted)) {
		stop("The two groups have different items omitted. Sum-score tests cannot be combined")
	}

	cb <- list(observed=e1$observed + e2$observed,
		   expected=e1$expected + e2$expected,
		   n=e1$n + e2$n,
		   omitted=e1$omitted)
	cb <- sumScoreEAPTestInternal(cb)
	cb
}

print.summary.sumScoreEAPTest <- function(x,...) {
	cat(sprintf("Latent distribution fit test (n=%d):\n", x$n))
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
sumScoreEAP <- function(grp, ..., qwidth=6.0, qpoints=49L, .twotier=TRUE) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}

	if (missing(qwidth) && !is.null(grp$qwidth)) { qwidth <- grp$qwidth }
	if (missing(qpoints) && !is.null(grp$qpoints)) { qpoints <- grp$qpoints }

	tbl <- ssEAP(grp, qwidth, qpoints, twotier=.twotier)
	rownames(tbl) <- 0:(nrow(tbl)-1)
	tbl
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
##' If you have missing data then you must specify
##' \code{minItemsPerScore}.  This option will set scores to NA when
##' there are too few items to make an accurate score estimate.  If
##' you are using the scores as point estimates without considering
##' the standard error then you should set \code{minItemsPerScore} as
##' high as you can tolerate. This will increase the amount of missing
##' data but scores will be more accurate. If you are carefully
##' considering the standard errors of the scores then you can set
##' \code{minItemsPerScore} to 1. This will mimic the behavior of most
##' other IFA software wherein scores are estimated if there is at
##' least 1 non-NA item for the score. However, it may make more sense
##' to set \code{minItemsPerScore} to 0. When set to 0, all NA rows
##' are scored to the prior distribution.
##'
##' @param grp a list with spec, param, and data
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param naAction deprecated, will be removed in the next release
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
EAPscores <- function(grp, ..., naAction=NULL) {
	if (length(list(...)) > 0) {
		stop(paste("Remaining parameters must be passed by name", deparse(list(...))))
	}

	if (!missing(naAction)) warning("naAction is deprecated")

	.Call(eap_wrapper, grp)
}

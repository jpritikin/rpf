##' Calculate cell moments
##' 
##' Popular moments include 2 (variance) and 4 (kurtosis).
##' 
##' @param spec list of item models
##' @param params data frame of item parameters, 1 per row
##' @param scores model derived person scores
##' @param m which moment
##' @return moment matrix
##' @docType methods
##' @export
rpf.1dim.moment <- function(spec, params, scores, m) {
  out <- array(dim=c(length(scores), length(spec)))
  for (ix in 1:length(spec)) {
    i <- spec[[ix]]
    prob <- rpf.prob(i, params[ix,], scores)
    Escore <- apply(prob, 1, function(r) sum(r * 0:(i@numOutcomes-1)))
    grid <- t(array(0:(i@numOutcomes-1), dim=c(i@numOutcomes, length(scores))))
    out[,ix] <- apply((grid - Escore)^m * prob, 1, sum)
  }
  out
}

##' Calculate residuals
##'
##' @param spec list of item models
##' @param params data frame of item parameters, 1 per row
##' @param responses persons in rows and items in columns
##' @param scores model derived person scores
##' @return residuals
##' @docType methods
##' @export
rpf.1dim.residual <- function(spec, params, responses, scores) {
  Zscore <- array(dim=c(length(scores), length(spec)))
  for (ix in 1:length(spec)) {
    i <- spec[[ix]]
    prob <- rpf.prob(i, params[ix,], scores)
    Escore <- apply(prob, 1, function(r) sum(r * 0:(i@numOutcomes-1)))
    data <- responses[,ix]
    if (!is.ordered(data)) { stop(paste("Column",ix,"is not an ordered factor")) }
    data <- unclass(data) - 1
    Zscore[,ix] <- data - Escore
  }
  Zscore
}

##' Calculate standardized residuals
##'
##' @param spec list of item models
##' @param params data frame of item parameters, 1 per row
##' @param responses persons in rows and items in columns
##' @param scores model derived person scores
##' @return standardized residuals
##' @docType methods
##' @export
rpf.1dim.stdresidual <- function(spec, params, responses, scores) {
  res <- rpf.1dim.residual(spec, params, responses, scores)
  variance <- rpf.1dim.moment(spec, params, scores, 2)
  res / sqrt(variance)
}

##' Calculate item and person fit statistics
##'
##' Exact distributional properties of these statistics are unknown
##' (Masters & Wright, 1997, p. 112).  For details on the calculation,
##' refer to Wright & Masters (1982, p. 100).
##'
##' The Wilson-Hilferty transformation is biased for less than 25 items.
##' To adjust Z scores for fewer items use wh.exact=FALSE.
##'
##' @param spec list of item models
##' @param params data frame of item parameters, 1 per row
##' @param responses persons in rows and items in columns
##' @param scores model derived person scores
##' @param margin for people 1, for items 2
##' @param na.rm remove NAs (default TRUE)
##' @param wh.exact whether to use the exact Wilson-Hilferty transformation (default TRUE)
##' @references Masters, G. N. & Wright, B. D. (1997). The Partial
##' Credit Model. In W. van der Linden & R. K. Kambleton (Eds.),
##' \emph{Handbook of modern item response theory}
##' (pp. 101-121). Springer.
##' 
##' Wilson, E. B., & Hilferty, M. M. (1931). The distribution of
##' chi-square. \emph{Proceedings of the National Academy of Sciences of the
##' United States of America,} 17, 684-688.
##' 
##' Wright, B. D. & Masters, G. N. (1982). \emph{Rating Scale
##' Analysis.} Chicago: Mesa Press.
##' @export
rpf.1dim.fit <- function(spec, params, responses, scores, margin, na.rm=TRUE, wh.exact=TRUE) {
# why permit na.rm=FALSE ? TODO
  r.var <- rpf.1dim.moment(spec, params, scores,2)
  r.k <- rpf.1dim.moment(spec, params, scores,4)
  r.z <- rpf.1dim.stdresidual(spec, params, responses, scores)
  wms.var <- apply(r.k - r.var^2, margin, sum, na.rm=na.rm)/
    apply(r.var, margin, sum, na.rm=na.rm)^2
  wms.sd <- sqrt(wms.var)
  fudge <- wms.sd/3
  if (!wh.exact) fudge <- 0
  outfit <- apply(r.z^2, margin, sum, na.rm=na.rm)/
                       apply(r.z, margin, function (l) sum(!is.na(l)))
  outfit.z <- (outfit^(1/3) - 1)*(3/wms.sd) + fudge
  infit <- apply(r.z^2 * r.var, margin, sum, na.rm=na.rm)/
                     apply(r.var, margin, sum, na.rm=na.rm)
  infit.z <- (infit^(1/3) - 1)*(3/wms.sd) + fudge
  df <- data.frame(infit, infit.z, outfit, outfit.z)
  if (margin == 2) {
    df$name <- rownames(params)    
  } else {
    df$name <- rownames(responses)
  }
  df
}

##' Find the point where an item provides mean maximum information
##'
##' @param spec an item spec
##' @param iparam an item parameter vector
##' @param grain the step size for numerical integration (optional)
rpf.mean.info1 <- function(spec, iparam, grain=.1) {
  range <- 9
  dim <- spec@dimensions
  if (dim != 1) stop("Not implemented")
  grid <- seq(-range, range, grain)
  info <- rpf.info(spec, iparam, grid)
  sum(info * grid) / sum(info)
}

##' Find the point where an item provides mean maximum information
##'
##' This is a point estimate of the mean difficulty of items that do
##' not offer easily interpretable parameters such as the Generalized
##' PCM. Since the information curve may not be unimodal, this
##' function integrates across the latent space.
##' 
##' @param spec list of item specs
##' @param param list or matrix of item parameters
##' @param grain the step size for numerical integration (optional)
rpf.mean.info <- function(spec, param, grain=.1) {
  ret <- list()
  for (ix in 1:length(spec)) {
    iparam <- c()
    if (is.list(param)) {
      iparam <- param[[ix]]
    } else {
      iparam <- param[ix,]
    }
    ret[[ix]] <- rpf.mean.info1(spec[[ix]], iparam, grain)
  }
  ret
}

##' Compute S-Chi-squared fit statistic for 1 item
##'
##' Implements the Kang & Chen (2007) polytomous extension to
##' S-Chi-squared statistic of Orlando & Thissen (2000). Fails in the
##' presence of missing data.
##'
##' NOTE: IRTPRO 2.1 uses an equal interval quadrature rule by default.
##'
##' WARNING: The algorithm for collapsing low-count cells has not been
##' tested thoroughly.
##'
##' @param spec a list of item specifications
##' @param param item paramters
##' @param free a matrix of the same shape as \code{param} indicating whether the parameter is free (TRUE) or fixed
##' @param item the item of interest
##' @param observed a matrix of observed raw scores by the outcome of the item of interest
##' @param quad the quadrature rule (default is G-H with 49 points)
##'
##' @references Kang, T. and Chen, T. T. (2007). An investigation of
##' the performance of the generalized S-Chisq item-fit index for
##' polytomous IRT models. ACT Research Report Series.
##'
##' Orlando, M. and Thissen, D. (2000). Likelihood-Based
##' Item-Fit Indices for Dichotomous Item Response Theory Models.
##' \emph{Applied Psychological Measurement, 24}(1), 50-64.
rpf.ot2000.chisq1 <- function(spec, param, free, item, observed, quad=NULL) {
  if (missing(quad)) quad <- rpf.GaussHermiteData(49)
  c.spec <- lapply(spec, function(m) {
    if (length(m@spec)==0) { stop("Not implemented") }
    else { m@spec }
  })

  out <- .Call(orlando_thissen_2000_wrapper, c.spec, param, item, observed, quad)
  observed <- c(out$observed)
  expected <- c(out$expected)
  mask <- expected!=0
  out$statistic <- sum((observed[mask] - expected[mask])^2 / expected[mask])
  out$df <- out$df - free;
  out$p.value <- dchisq(out$statistic, out$df)
  out
}

##' Compute S-Chi-squared fit statistic for a set of items
##'
##' Runs \code{\link{rpf.ot2000.chisq1}} for every item and accumulates
##' the results.
##'
##' @param spec a list of item specifications
##' @param param item paramters in columns
##' @param free a matrix of the same shape as \code{param} indicating whether the parameter is free (TRUE) or fixed
##' @param data a data frame or matrix of response patterns, 1 per row
##' @param quad the quadrature rule (default is G-H with 49 points)
##' @return
##' a list of output from \code{\link{rpf.ot2000.chisq1}}
rpf.ot2000.chisq <- function(spec, param, free, data, quad=NULL) {
  if (missing(quad)) quad <- rpf.GaussHermiteData(49)
  if (is.data.frame(data)) {
    data <- sapply(data, unclass)
  }
  if (dim(data)[2] != length(spec)) stop("Dim mismatch between data and spec")
  if (dim(param)[2] != length(spec)) stop("Dim mismatch between param and spec")
  context <- 1:length(spec)
  got = list()
  for (interest in 1:length(spec)) {
    if (0) {
      # hard to persuade R to do this correctly
      ob.table <- table(apply(data[,context[context != interest]], 1, sum), data[,interest])
      sumscores <- (length(spec)-1):(sum(sapply(spec, function(m) slot(m,'numOutcomes'))) - spec[[interest]]@numOutcomes)
      observed <- array(0, dim=c(length(sumscores), spec[[interest]]@numOutcomes))
      rowmap <- match(sumscores, rownames(ob.table))
      rowmap <- rowmap[!is.na(rowmap)]
      observed[rowmap,] <- ob.table
    } else {
      observed <- .Call(sumscore_observed, sum(sapply(spec, function(m) slot(m,'numOutcomes'))),
                        data, interest, spec[[interest]]@numOutcomes)
    }
    ot.out <- rpf.ot2000.chisq1(spec, param, sum(free[,interest]), interest, observed, quad)
    got[[interest]] <- ot.out
  }
  got
}

##' Compute Gauss-Hermite quadrature rule
##'
##' Computes Gauss-Hermite quadrature rule of requested order using the
##' Golub-Welsch algorithm. This is very fast and numerically
##' stable. It can handle quadrature of order 1000+.
##' @param
##' n Order of Gauss-Hermite rule to compute (number of nodes)
##' @return A list containing the node positions (x) and the
##' quadrature weights (w) for the requested rule.
##' @author
##' Alexander W Blocker <ablocker@@gmail.com>
##' @aliases
##' rpf_GaussHermiteData
##' @references
##' Golub, G. H. and Welsch, J. H. (1969). Calculation of
##' Gauss Quadrature Rules. Mathematics of Computation 23 (106):
##' 221-230
rpf.GaussHermiteData <- function(n) .Call(rpf_GaussHermiteData, n)

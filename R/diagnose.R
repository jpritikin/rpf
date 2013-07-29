##' Calculate cell central moments
##' 
##' Popular central moments include 2 (variance) and 4 (kurtosis).
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
    prob <- t(rpf.prob(i, params[,ix], scores))  # remove t() TODO
    Escore <- apply(prob, 1, function(r) sum(r * 0:(i@outcomes-1)))
    grid <- t(array(0:(i@outcomes-1), dim=c(i@outcomes, length(scores))))
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
    prob <- t(rpf.prob(i, params[,ix], scores))  # remove t() TODO
    Escore <- apply(prob, 1, function(r) sum(r * 0:(i@outcomes-1)))
    data <- responses[,ix]
    if (is.ordered(data)) {
      data <- unclass(data) - 1
    } else if (is.factor(data)) {
      stop(paste("Column",ix,"is an unordered factor"))
    }
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

##' Calculate item and person Rasch fit statistics
##'
##' Note: These statistics are only appropriate if all discrimination
##' parameters are fixed equal and there is no missing data.
##'
##' Exact distributional properties of these statistics are unknown
##' (Masters & Wright, 1997, p. 112).  For details on the calculation,
##' refer to Wright & Masters (1982, p. 100).
##'
##' The Wilson-Hilferty transformation is biased for less than 25 items.
##' Consider wh.exact=FALSE for less than 25 items.
##'
##' @param spec list of item models
##' @param params matrix of item parameters, 1 per column
##' @param responses persons in rows and items in columns
##' @param scores model derived person scores
##' @param margin for people 1, for items 2
##' @param wh.exact whether to use the exact Wilson-Hilferty transformation
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
rpf.1dim.fit <- function(spec, params, responses, scores, margin, wh.exact=TRUE) {
  if (any(is.na(responses))) warning("Rasch fit statistics should not be used with missing data")  # true? TODO

  if (dim(params)[2] < 25 && wh.exact) {
    if (missing(wh.exact)) {
      wh.exact <- FALSE
      warning("Consider wh.exact=FALSE for less than 25 items")
    }
  }
  if (dim(params)[2] > 25 && !wh.exact) {
    if (missing(wh.exact)) {
      wh.exact <- TRUE
      warning("Consider wh.exact=TRUE for more than 25 items")
    }
  }

  exclude.col <- c()
  outcomes <- sapply(spec, function(s) s@outcomes)
  for (ix in 1:dim(responses)[2]) {
    kat <- sum(table(responses[,ix]) > 0)
    if (kat != outcomes[ix]) {
      exclude.col <- c(exclude.col, ix)
      warning(paste("Excluding item", colnames(responses)[ix], "because outcomes !=", outcomes[ix]))
    }
  }

  if (length(exclude.col)) {
    responses <- responses[,-exclude.col]
    spec <- spec[-exclude.col]
    params <- params[,-exclude.col]
    outcomes <- outcomes[-exclude.col]
  }

  exclude.row <- c()
  for (ix in 1:dim(responses)[1]) {
    r1 <- sapply(responses[ix,], unclass)
    if (any(is.na(r1))) next
    if (all(r1 == 1) || all(r1 == outcomes)) {
      exclude.row <- c(exclude.row, ix)
      warning(paste("Excluding response", rownames(responses)[ix], "because minimum or maximum"))
    }
  }

  if (length(exclude.row)) {
    responses <- responses[-exclude.row,]
    scores <- scores[-exclude.row]
  }

  na.rm=TRUE
  r.z <- rpf.1dim.stdresidual(spec, params, responses, scores)
  r.var <- rpf.1dim.moment(spec, params, scores,2)
  r.var[is.na(r.z)] <- NA
  r.k <- rpf.1dim.moment(spec, params, scores,4)
  r.k[is.na(r.z)] <- NA

  outfit.var <- r.var
  outfit.var[r.var^2 < 1e-5] <- sqrt(1e-5)
  outfit.n <- apply(r.var, margin, function(l) sum(!is.na(l)))
  outfit.sd <- sqrt(apply(r.k / outfit.var^2, margin, sum, na.rm=na.rm) /
                    outfit.n^2 - 1/outfit.n)
  outfit.sd[outfit.sd > 1.4142] <- 1.4142
  outfit.fudge <- outfit.sd/3

  infit.sd <- sqrt(apply(r.k - r.var^2, margin, sum, na.rm=na.rm)/
    apply(r.var, margin, sum, na.rm=na.rm)^2)
  infit.sd[infit.sd > 1.4142] <- 1.4142
  infit.fudge <- infit.sd/3
  if (!wh.exact) {
    infit.fudge <- 0
    outfit.fudge <- 0
  }

  outfit <- apply(r.z^2, margin, sum, na.rm=na.rm)/
                       apply(r.z, margin, function (l) sum(!is.na(l)))
  outfit.z <- (outfit^(1/3) - 1)*(3/outfit.sd) + outfit.fudge

  infit <- apply(r.z^2 * r.var, margin, sum, na.rm=na.rm)/
                     apply(r.var, margin, sum, na.rm=na.rm)
  infit.z <- (infit^(1/3) - 1)*(3/infit.sd) + infit.fudge

  df <- data.frame(n=outfit.n, infit, infit.z, outfit, outfit.z)
  if (margin == 2) {
    df$name <- colnames(params)
  } else {
    df$name <- rownames(responses)
  }
  df
}

##' Find the point where an item provides mean maximum information
##'
##' WARNING: This function is experimental and may disappear.
##'
##' @param spec an item spec
##' @param iparam an item parameter vector
##' @param grain the step size for numerical integration (optional)
rpf.mean.info1 <- function(spec, iparam, grain=.1) {
  range <- 9
  dim <- spec@factors
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
##' @param quad the quadrature rule (default is equally spaced intervals with 49 points)
##'
##' @references Kang, T. and Chen, T. T. (2007). An investigation of
##' the performance of the generalized S-Chisq item-fit index for
##' polytomous IRT models. ACT Research Report Series.
##'
##' Orlando, M. and Thissen, D. (2000). Likelihood-Based
##' Item-Fit Indices for Dichotomous Item Response Theory Models.
##' \emph{Applied Psychological Measurement, 24}(1), 50-64.
rpf.ot2000.chisq1 <- function(spec, param, free, item, observed, quad=NULL) {
  if (missing(quad)) {
    n <- 49
    width <- 6
    x <- seq(-width, width, length.out=n)
    quad <- list(x=x, w=dnorm(x)/sum(dnorm(x)))
  }
  c.spec <- lapply(spec, function(m) {
    if (length(m@spec)==0) { stop("Item model",m,"is not fully implemented") }
    else { m@spec }
  })

  max.param <- max(vapply(spec, rpf.numParam, 0))
  if (dim(param)[1] < max.param) {
    stop(paste("param matrix must have", max.param ,"rows"))
  }
  out <- .Call(orlando_thissen_2000_wrapper, c.spec, param, item, observed, quad)
  out$orig.observed <- out$observed
  out$orig.expected <- out$expected
  kc <- .Call(kang_chen_2007_wrapper, out$orig.observed, out$orig.expected)
  out$observed <- observed <- kc$observed
  out$expected <- expected <- kc$expected
  mask <- expected!=0
  out$statistic <- sum((observed[mask] - expected[mask])^2 / expected[mask])
  out$df <- out$df - free - kc$collapsed;
  out$p.value <- dchisq(out$statistic, out$df)
  out
}

##' Compute S-Chi-squared fit statistic for a set of items
##'
##' Runs \code{\link{rpf.ot2000.chisq1}} for every item and accumulates
##' the results.
##'
##' TODO: Handle multiple factors with mean and covariance parameters.
##'
##' @param spec a list of item specifications
##' @param param item paramters in columns
##' @param free a matrix of the same shape as \code{param} indicating whether the parameter is free (TRUE) or fixed
##' @param data a data frame or matrix of response patterns, 1 per row
##' @param quad the quadrature rule (default is equally spaced intervals with 49 points)
##' @return
##' a list of output from \code{\link{rpf.ot2000.chisq1}}
rpf.ot2000.chisq <- function(spec, param, free, data, quad=NULL) {
  if (missing(quad)) {
    n <- 49
    width <- 6
    x <- seq(-width, width, length.out=n)
    quad <- list(x=x, w=dnorm(x)/sum(dnorm(x)))
  }
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
      sumscores <- (length(spec)-1):(sum(sapply(spec, function(m) slot(m,'outcomes'))) - spec[[interest]]@outcomes)
      observed <- array(0, dim=c(length(sumscores), spec[[interest]]@outcomes))
      rowmap <- match(sumscores, rownames(ob.table))
      rowmap <- rowmap[!is.na(rowmap)]
      observed[rowmap,] <- ob.table
    } else {
      observed <- .Call(sumscore_observed, sum(sapply(spec, function(m) slot(m,'outcomes'))),
                        data, interest, spec[[interest]]@outcomes)
    }
    ot.out <- rpf.ot2000.chisq1(spec, param, sum(free[,interest]), interest, observed, quad)
    got[[interest]] <- ot.out
  }
  got
}

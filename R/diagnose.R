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
    if (length(data) != length(Escore)) stop("Length mismatch")
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
##' parameters are fixed equal and items are conditionally independent
##' (see \code{\link{chen.thissen.1997}}).  A best effort is made to
##' cope with missing data.
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
##' @param group spec, params, data, and scores can be provided in a list instead of as arguments
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
rpf.1dim.fit <- function(spec, params, responses, scores, margin, group=NULL, wh.exact=TRUE) {
    if (!missing(group)) {
        spec <- group$spec
        params <- group$param
        responses <- group$data
        scores <- group$scores[,1]  # should not assume first score TODO
    }

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

thetaComb <- function(theta, nfact)  #copied from mirt
{
  if (nfact == 1L){
    Theta <- matrix(theta)
  } else {
    thetalist <- vector('list', nfact)
    for(i in 1L:nfact)
      thetalist[[i]] <- theta
    Theta <- as.matrix(expand.grid(thetalist))
  }	
  return(Theta)
}

##' Compute the ordinal gamma association statistic
##'
##' @param mat a cross tabulation matrix
##' @references
##' Agresti, A. (1990). Categorical data analysis. New York: Wiley.
ordinal.gamma <- function(mat) .Call(ordinal_gamma_wrapper, mat)

# root mean squared statistic (sqrt omitted)
ms <- function(observed, expected, draws) {
  draws * sum((observed - expected)^2)
}

P.cdf.fn <- function(x, g.var, t) {
  sapply(t, function (t1) {
    n <- length(g.var)
    num <- exp(1-t1) * exp(1i * t1 * sqrt(n))
    den <- pi * (t1 - 1/(1-1i*sqrt(n)))
    pterm <- prod(sqrt(1 - 2*(t1-1)*g.var/x + 2i*t1*g.var*sqrt(n)/x))
    Im(num / (den * pterm))
  })
}

##' Compute the P value that the observed and expected tables come from the same distribution
##'
##' This method dramatically improves upon Pearson's X^2
##' goodness-of-fit test.  In contrast to Pearson's X^2, no ad hoc cell
##' collapsing is needed to avoid an inflated false positive rate.
##' The statistic rapidly converges to the Monte-Carlo estimate
##' as the number of draws increases. In contrast to Pearson's
##' X^2, the order of the matrices doesn't matter. This test is
##' commutative with respect to its arguments.
##' 
##' @param observed observed matrix
##' @param expected expected matrix
##' @return The P value indicating whether the two tables come from
##' the same distribution. For example, a significant result (P <
##' alpha level) rejects the hypothesis that the two matrices are from
##' the same distribution.
##' @references Perkins, W., Tygert, M., & Ward, R. (2011). Computing
##' the confidence levels for a root-mean-square test of
##' goodness-of-fit. \emph{Applied Mathematics and Computations,
##' 217}(22), 9072-9084.
##' @examples
##' draws <- 17
##' observed <- matrix(c(.294, .176, .118, .411), nrow=2) * draws
##' expected <- matrix(c(.235, .235, .176, .353), nrow=2) * draws
##' ptw2011.gof.test(observed, expected)  # not signficiant

ptw2011.gof.test <- function(observed, expected) {
  orig.draws <- sum(observed)
  observed <- observed / orig.draws
  expected <- expected / orig.draws

  X <- ms(observed, expected, orig.draws)
  if (X == 0) return(1)

  n <- prod(dim(observed))
  D <- diag(1/c(expected))
  P <- matrix(-1/n, n,n)
  diag(P) <- 1 - 1/n
  B <- P %*% D %*% P
  g.var <- 1 / eigen(B, only.values=TRUE)$values[-n]
  
  # Eqn 8 needs n variances, but matrix B (Eqn 6) only has n-1 non-zero
  # eigenvalues. Perhaps n-1 degrees of freedom?
  
#  plot(function(t) P.cdf.fn(X, g.var, t), 0, 40)
  # 310 points should be good enough for 500 bins
  # Perkins, Tygert, Ward (2011, p. 10)
  got <- integrate(function(t) P.cdf.fn(X, g.var, t), 0, 40, subdivisions=310L)
  p.value <- 1 - got$value
  smallest <- 6.3e-16  # approx exp(-35)
  if (p.value < smallest) p.value <- smallest
  p.value
}

##' Computes local dependence indices for all pairs of items
##'
##' Item Factor Analysis makes two assumptions: (1) that the latent
##' distribution is reasonably approximated by the multivariate Normal
##' and (2) that items are conditionally independent. This test
##' examines the second assumption. The presence of locally dependent
##' items can inflate the precision of estimates causing a test to
##' seem more accurate than it really is.
##'
##' Statically significant entries suggest that the item pair has
##' local dependence. Since log(.01)=-4.6, an absolute magitude of 5
##' is a reasonable cut-off. Positive entries indicate that the two
##' items are more correlated than expected. These items may share an
##' unaccounted for latent dimension. Consider a redesign of the items
##' or the use of testlets for scoring. Negative entries indicate that
##' the two items are less correlated than expected.
##'
##' @param grp a list with the spec, param, mean, and cov describing the group
##' @param data data
##' @param inames a subset of items to examine
##' @param qwidth quadrature width
##' @param qpoints number of equally spaced quadrature points
##' @param method method to use to calculate P values. The default
##' ("rms") uses the root mean square statistic (see \code{\link{ptw2011.gof.test}}).
##' To obtain the traditional Pearson X^2 statistic, use method="pearson".
##' @return a lower triangular matrix of log P values with the sign
##' determined by relative association between the observed and
##' expected tables (see \code{\link{ordinal.gamma}})
##' @references Chen, W.-H. & Thissen, D. (1997). Local dependence
##' indexes for item pairs using Item Response Theory. \emph{Journal
##' of Educational and Behavioral Statistics, 22}(3), 265-289.
##'
##' Wainer, H. & Kiely, G. L. (1987). Item clusters and computerized
##' adaptive testing: A case for testlets.  \emph{Journal of
##' Educational measurement, 24}(3), 185--201.
chen.thissen.1997 <- function(grp, data=NULL, inames=NULL, qwidth=6, qpoints=49, method="rms") {
  if (is.null(colnames(grp$param))) stop("Item parameter columns must be named")

  if (missing(data)) {
      data <- grp$data
  }
  if (method != "rms" && method != "pearson") stop(paste("Unknown method", method))
  if (missing(inames)) {
    inames <- colnames(grp$param)
  }
  if (length(inames) < 2) stop("At least 2 items are required")

  spec <- grp$spec
  if (length(spec) < dim(grp$param)[2]) {
      if (dim(grp$param)[2] %% length(spec) != 0) stop("Length of spec must match # of items")
      rep <- dim(grp$param)[2] %/% length(spec)
      while (rep > 1) {
          spec <- c(spec, grp$spec)
          rep <- rep - 1
      }
  }

  if (!is.data.frame(data)) {
      data <- as.data.frame(data)  #safe? TODO
  }

  # assume param and data in the same order? TODO
  items <- match(inames, colnames(grp$param))

  result <- list()
  gamma <- matrix(NA, length(items), length(items))
  dimnames(gamma) <- list(inames, inames)
  raw <- matrix(NA, length(items), length(items))
  dimnames(raw) <- list(inames, inames)
  std <- matrix(NA, length(items), length(items))
  dimnames(std) <- list(inames, inames)
  pval <- matrix(NA, length(items), length(items))
  dimnames(pval) <- list(inames, inames)

  # sort item pairs by fset TODO
  for (iter1 in 2:length(items)) {
    for (iter2 in 1:(iter1-1)) {
      i1 <- items[iter1]
      i2 <- items[iter2]
      observed <- table(data[,c(i1,i2)])
      N <- sum(observed)
      s1 <- spec[[i1]]
      s2 <- spec[[i2]]
      expected <- matrix(NA, s1@outcomes, s2@outcomes)
      dimnames(expected) <- dimnames(observed)

      # only integrate over relevant factors
      f1 <- which(grp$param[1:s1@factors,i1] != 0)
      f2 <- which(grp$param[1:s2@factors,i2] != 0)
      fset <- union(f1, f2)
      mean <- grp$mean[fset]
      cov <- grp$cov[fset,fset, drop=FALSE]
      theta <- thetaComb(seq(-qwidth,qwidth,length.out=qpoints), length(mean))
      prior <- mvtnorm::dmvnorm(theta, mean, cov)
      prior <- prior/sum(prior)
      palette <- t(theta)
      th1 <- palette[match(f1, fset),, drop=FALSE]
      th2 <- palette[match(f2, fset),, drop=FALSE]

      p1 <- rpf.prob(s1, grp$param[,i1], th1)
      p2 <- rpf.prob(s2, grp$param[,i2], th2)
      for (o1 in 1:s1@outcomes) {
        for (o2 in 1:s2@outcomes) {
          expected[o1,o2] <- N * sum(p1[o1,] * p2[o2,] * prior)
        }
      }
      s <- ordinal.gamma(observed) - ordinal.gamma(expected)
      if (!is.finite(s) || is.na(s) || s==0) s <- 1
      info <- list(observed=observed, expected=expected, sign=sign(s), gamma=s)
      gamma[iter1, iter2] <- s

      if (method == "rms") {
          raw[iter1, iter2] <- ms(observed, expected, sum(observed))
          pval[iter1, iter2] <- sign(s) * -log(ptw2011.gof.test(observed, expected))
      } else {
          x2 <- sum((observed - expected)^2 / expected)
          df <- (s1@outcomes-1) * (s2@outcomes-1)
          info <- c(info, x2=x2, df=df)

          raw[iter1, iter2] <- sign(s) * x2
          std[iter1, iter2] <- sign(s) * abs((x2 - df)/sqrt(2*df))
          pval[iter1, iter2] <- sign(s) * -pchisq(x2, df, lower.tail=FALSE, log.p=TRUE)
      }

      result[[paste(inames[iter1], inames[iter2], sep=":")]] <- info
    }
  }
  list(pval=pval, std=std, raw=raw, gamma=gamma, detail=result)
}


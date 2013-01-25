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

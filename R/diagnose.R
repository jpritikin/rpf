##' Calculate standardized residuals
##'
##' Squaring and averaging the margins produce item and person outfit
##' statistics. For details, see Embretson & Reise (2000, pp. 237-238).
##'
##' @param items list of item models
##' @param params data frame of item parameters, 1 per row
##' @param responses persons in rows and items in columns
##' @param scores model derived person scores
##' @return standardized residuals
##' @docType methods
##' @aliases
##' rpf.residuals,rpf.1dim,list,data.frame,data.frame,numeric-method
##' @export
##' @references Embretson, S. E. & Reise S. P. (2000) Item response
##' theory for psychologists. Lawrence Erlbaum.
rpf.1dim.residuals <- function(items, params, responses, scores) {
  Zscore <- array(dim=c(length(scores), length(items)))
  for (ix in 1:length(items)) {
    i <- items[[ix]]
    espt.prob <- (rpf.prob(i, params[ix,], scores))
    Escore <- apply(espt.prob, 1, function(r) sum(r * 1:i@numOutcomes))
    Vscore <- numeric(length(scores))
    data <- responses[,ix]
    if (!is.ordered(data)) { stop(paste("Column",ix,"is not an ordered factor")) }
    data <- unclass(data)
    for (sx in 1:length(scores)) {
      Vscore[sx] <- sum((1:i@numOutcomes - data[sx])^2 * espt.prob[sx,])
    }
    Zscore[,ix] <- (data - Escore) / sqrt(Vscore)
  }
  Zscore
}

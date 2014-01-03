##' Unpack a two-tier model
##'
##' A design matrix can be used to compress a bifactor or two-tier
##' model into primary factors and a partition of items into specific
##' factors. This function unpacks the design matrix such that the
##' specific factors are given their own separate loading across all
##' the items.
##' 
##' @param grp a list containing the design matrix, spec, param, mean, cov, scores, and data
##' @return the same model without a design matrix
##' @examples
##' # Need example here, TODO
unpack.2tier <- function(grp) {
  nfact <- max(grp$design, na.rm=TRUE)
  if (nfact <= nrow(grp$design)) return(grp)
  param <- matrix(0, nrow=nfact - nrow(grp$design) + nrow(grp$param), ncol=ncol(grp$param))
  param[(nfact+1):nrow(param),] <- grp$param[(nrow(grp$design)+1):nrow(grp$param),]
  colnames(param) <- colnames(grp$param)
  for (d in 1:nfact) {
    mask <- grp$design==d
    mask <- mask & !is.na(mask)
    param[d, apply(mask, 2, any)] <- grp$param[1:nrow(grp$design),][mask]
  }
  spec <- sapply(grp$spec, rpf.modify, nfact)
  names(spec) <- names(grp$spec)
  list(spec=spec,
       param=param,
       mean=grp$mean,
       cov=grp$cov,
       scores=grp$scores,
       data=grp$data)
}

sumScoreEAP <- function(grp, width, pts) {
	if (missing(width)) width <- 6
	if (missing(pts)) pts <- 49L
	.Call(ssEAP_wrapper, grp, width, pts)
}

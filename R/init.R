.onLoad <- function(libname, pkgname) {
	cores <- detectCores(logical=TRUE)
	.Call(setNumberOfCores, cores)
	options("mc.cores"=cores)  # maybe the wrong way to set the default? TODO
}

.onAttach <- function(libname, pkgname) {
	if (! .Call(hasOpenMP_wrapper)) {
		packageStartupMessage("RPF is not compiled to take advantage of computers with multiple cores.")
	}
}

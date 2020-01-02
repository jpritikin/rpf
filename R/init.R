.onLoad <- function(libname, pkgname) {
	if ("package:parallel" %in% search()) {
		cores <- detectCores(logical=TRUE)
		.Call(`_rpf_setNumberOfCores`, cores)
		#packageStartupMessage(paste("OpenMP will use", cores, "cores"))
	}
}

.onAttach <- function(libname, pkgname) {
	if (! .Call(`_rpf_has_openmp`)) {
		packageStartupMessage("RPF is not compiled to take advantage of computers with multiple cores.")
	}
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call('_rpf_registerCCallable', PACKAGE = 'rpf')
})

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "../inst/include/libirt-rpf.h"

void getMatrixDims(SEXP r_theta, int *rows, int *cols);
SEXP orlando_thissen_2000(SEXP r_spec, SEXP r_param, SEXP r_item, SEXP r_observed, SEXP r_quad);
SEXP omxGaussHermiteData(SEXP points);
SEXP sumscore_observed(SEXP r_high, SEXP r_data, SEXP r_interest, SEXP r_outcomes);

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "../inst/include/libifa-rpf.h"

void getMatrixDims(SEXP r_theta, int *rows, int *cols);
SEXP orlando_thissen_2000(SEXP r_spec, SEXP r_param, SEXP r_item, SEXP r_observed, SEXP r_quad);
SEXP sumscore_observed(SEXP r_high, SEXP r_data, SEXP r_interest, SEXP r_outcomes);
SEXP kang_chen_2007_wrapper(SEXP r_observed_orig, SEXP r_expected_orig);
SEXP gamma_cor(SEXP r_mat);

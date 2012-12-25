#include <R.h>
#include <Rinternals.h>
#include "libirt-rpf.h"

SEXP rpf_1dim_drm_logprob_wrapper(SEXP r_param, SEXP r_theta)
{
  if (length(r_param) != 3) error("Wrong parameter length");

  double *theta = REAL(r_theta);
  int thetaLen = length(r_theta);

  const int numOutcomes = 2;
  SEXP outsxp;
  PROTECT(outsxp = allocMatrix(REALSXP, numOutcomes, thetaLen));
  double *out = REAL(outsxp);

  for (int px=0; px < thetaLen; px++) {
    double th = theta[px];
    irt_rpf_1dim_drm_logprob(REAL(r_param), th, out+px*numOutcomes);
  }

  UNPROTECT(1);
  return outsxp;
}

static void
getMatrixDims(SEXP r_theta, int *rows, int *cols)
{
    SEXP matrixDims;
    PROTECT(matrixDims = getAttrib(r_theta, R_DimSymbol));
    int *dimList = INTEGER(matrixDims);
    *rows = dimList[0];
    *cols = dimList[1];
    UNPROTECT(1);
}

SEXP rpf_mdim_drm_logprob_wrapper(SEXP r_numDims,
				  SEXP r_param, SEXP r_theta)
{
  int numDims = asInteger(r_numDims);

  if (length(r_param) != numDims+2) error("Wrong parameter length");

  double *theta = REAL(r_theta);
  int numPersons;
  int numAbilities;
  getMatrixDims(r_theta, &numAbilities, &numPersons);
  if (numAbilities != numDims) error("Ability dimension mismatch");

  const int numOutcomes = 2;
  SEXP outsxp;
  PROTECT(outsxp = allocMatrix(REALSXP, numOutcomes, numPersons));
  double *out = REAL(outsxp);

  for (int px=0; px < numPersons; px++) {
    irt_rpf_mdim_drm_logprob(numDims, REAL(r_param), theta + px*numAbilities,
			     out+px*numOutcomes);
  }

  UNPROTECT(1);
  return outsxp;
}

SEXP rpf_1dim_gpcm_logprob_wrapper(SEXP r_numOutcomes,
				   SEXP r_param, SEXP r_theta)
{
  int numOutcomes = asInteger(r_numOutcomes);

  if (length(r_param) != numOutcomes) error("Wrong parameter length");

  double *theta = REAL(r_theta);
  int thetaLen = length(r_theta);

  SEXP outsxp;
  PROTECT(outsxp = allocMatrix(REALSXP, numOutcomes, thetaLen));
  double *out = REAL(outsxp);

  for (int px=0; px < thetaLen; px++) {
    double th = theta[px];
    irt_rpf_1dim_gpcm_logprob(numOutcomes, REAL(r_param),
			      th, out+px*numOutcomes);
  }

  UNPROTECT(1);
  return outsxp;
}

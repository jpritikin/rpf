#include <R.h>
#include <Rinternals.h>

static void
rpf_1dim_gpcm_prob(int numOutcomes, double D, double *param,
		   int thetaLen, double *theta, double *out)
{
  double discr = param[0] * D;
  double term1[numOutcomes];

  term1[numOutcomes - 1] = 0;

  for (int px=0; px < thetaLen; px++) {
    double th = theta[px];
    for (int tx=0; tx < numOutcomes-1; tx++) {
      term1[tx] = -discr * (th - param[tx+1]);
    }
    double sum = 0;
    double denom = 0;
    double term2[numOutcomes];
    for (int tx=numOutcomes-1; tx >= 0; tx--) {
      sum += term1[tx];
      double tmp = exp(sum);
      term2[tx] = tmp;
      denom += tmp;
    }
    for (int tx=0; tx < numOutcomes; tx++) {
      out[tx * thetaLen + px] = term2[tx] / denom;
    }
  }
}

SEXP rpf_1dim_gpcm_prob_wrapper(SEXP r_numOutcomes, SEXP r_D,
				SEXP r_param, SEXP r_theta)
{
  int numOutcomes = asInteger(r_numOutcomes);
  double D = asReal(r_D);

  if (length(r_param) != numOutcomes) {
    error("Expecting parameter vector of length %d, got %d",
	  numOutcomes, length(r_param));
  }

  double *theta = REAL(r_theta);
  if (isMatrix(r_theta)) {
    SEXP matrixDims;
    PROTECT(matrixDims = getAttrib(r_theta, R_DimSymbol));
    int *dimList = INTEGER(matrixDims);
    if (dimList[1] != 1) {
      error("Expecting 1 theta parameter per person, got %d",
	    dimList[1]);
    }
    UNPROTECT(1);
  } else if (isVector(r_theta)) {
    // ok
  } else {
    error("Expecting a matrix or vector");
  }

  int thetaLen = length(r_theta);

  SEXP out;
  PROTECT(out = allocMatrix(REALSXP, thetaLen, numOutcomes));

  rpf_1dim_gpcm_prob(numOutcomes, D, REAL(r_param),
		     thetaLen, theta, REAL(out));

  UNPROTECT(1);
  return out;
}

#include "rpf.h"

static SEXP
get_model_names(SEXP name)
{
  const char *target = CHAR(STRING_ELT(name, 0));
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, 1));
  REAL(ret)[0] = NA_REAL;
  for (int sx=0; sx < librpf_numModels; sx++) {
    if (strcmp(librpf_model[sx].name, target) == 0) {
      REAL(ret)[0] = sx;
    }
  }
  UNPROTECT(1);
  return ret;
}

void getMatrixDims(SEXP r_theta, int *rows, int *cols)
{
    SEXP matrixDims;
    PROTECT(matrixDims = getAttrib(r_theta, R_DimSymbol));
    int *dimList = INTEGER(matrixDims);
    *rows = dimList[0];
    *cols = dimList[1];
    UNPROTECT(1);
}

static SEXP
rpf_numSpec_wrapper(SEXP r_spec)
{
  if (length(r_spec) < RPF_ISpecCount)
    error("Item spec must be of length %d, not %d", RPF_ISpecCount, length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= librpf_numModels)
    error("Item model %d out of range", id);

  int numSpec = (*librpf_model[id].numSpec)(spec);

  SEXP ret;
  PROTECT(ret = allocVector(INTSXP, 1));
  INTEGER(ret)[0] = numSpec;
  UNPROTECT(1);

  return ret;
}

static SEXP
rpf_numParam_wrapper(SEXP r_spec)
{
  if (length(r_spec) < RPF_ISpecCount)
    error("Item spec must be of length %d, not %d", RPF_ISpecCount, length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= librpf_numModels)
    error("Item model %d out of range", id);

  int numParam = (*librpf_model[id].numParam)(spec);

  SEXP ret;
  PROTECT(ret = allocVector(INTSXP, 1));
  INTEGER(ret)[0] = numParam;
  UNPROTECT(1);

  return ret;
}

static SEXP
rpf_paramInfo_wrapper(SEXP r_spec, SEXP r_paramNum)
{
  if (length(r_spec) < RPF_ISpecCount)
    error("Item spec must be of length %d, not %d", RPF_ISpecCount, length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= librpf_numModels)
    error("Item model %d out of range", id);

  int pnum = asInteger(r_paramNum);
  int numParam = (*librpf_model[id].numParam)(spec);
  if (pnum < 0 || pnum >= numParam) error("Item model %d has %d parameters", id, numParam);

  int type;
  double upper, lower;
  (*librpf_model[id].paramInfo)(spec, pnum, &type, &upper, &lower);

  int len = 3;
  SEXP names, ans;
  PROTECT(names = allocVector(STRSXP, len));
  PROTECT(ans = allocVector(VECSXP, len));
  int lx = 0;
  SET_STRING_ELT(names, lx, mkChar("type"));
  SET_VECTOR_ELT(ans,   lx, ScalarInteger(type));
  SET_STRING_ELT(names, ++lx, mkChar("upper"));
  SET_VECTOR_ELT(ans,   lx, ScalarReal(isfinite(upper)? upper : NA_REAL));
  SET_STRING_ELT(names, ++lx, mkChar("lower"));
  SET_VECTOR_ELT(ans,   lx, ScalarReal(isfinite(lower)? lower : NA_REAL));
  namesgets(ans, names);
  UNPROTECT(2);

  return ans;
}

int unpack_theta(int dims, double *param, int numAbilities, double *theta, double *out)
{
  if (numAbilities == dims) {
    for (int dx=0; dx < dims; ++dx) {
      double th = theta[dx];
      if (!isfinite(th)) return 0;
      out[dx] = th;
    }
  } else {
    int ax = 0;
    for (int dx=0; dx < dims; ++dx) {
      if (param[dx] == 0) continue;
      double th = theta[ax]; // could read uninitialized memory, but we detect below
      if (!isfinite(th)) return 0;
      out[dx] = th;
      ++ax;
    }
    if (ax != numAbilities) {
      error("Item has %d nonzero dims but given %d abilities", ax, numAbilities);
    }
  }
  return 1;
}

static SEXP
rpf_prob_wrapper(SEXP r_spec, SEXP r_param, SEXP r_theta)
{
  if (length(r_spec) < RPF_ISpecCount)
    error("Item spec must be of length %d, not %d", RPF_ISpecCount, length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= librpf_numModels)
    error("Item model %d out of range", id);

  int numSpec = (*librpf_model[id].numSpec)(spec);
  if (length(r_spec) < numSpec)
    error("Item spec must be of length %d, not %d", numSpec, length(r_spec));
    
  int numParam = (*librpf_model[id].numParam)(spec);
  if (length(r_param) < numParam)
    error("Item has %d parameters, only %d given", numParam, length(r_param));

  const int numOutcomes = spec[RPF_ISpecOutcomes];
  const int dims = spec[RPF_ISpecDims];

  int numPeople;
  int numAbilities = 1;
  if (dims == 1) {
    numPeople = length(r_theta);
  } else {
    getMatrixDims(r_theta, &numAbilities, &numPeople);
  }

  SEXP outsxp;
  PROTECT(outsxp = allocMatrix(REALSXP, numOutcomes, numPeople));
  double *out = REAL(outsxp);
  double *theta = REAL(r_theta);
  double *param = REAL(r_param);
    
  for (int px=0; px < numPeople; px++) {
    double thBuf[dims];
    if (!unpack_theta(dims, param, numAbilities, theta + px*numAbilities, thBuf)) {
	for (int ox=0; ox < numOutcomes; ox++) {
	  out[px*numOutcomes + ox] = NA_REAL;
	}
	continue;
    }
    (*librpf_model[id].prob)(spec, param, thBuf, out+px*numOutcomes);
    for (int ox=0; ox < numOutcomes; ox++) {
      double prob = out[px*numOutcomes + ox];
      if (!isfinite(prob)) {
	out[px*numOutcomes + ox] = NA_REAL;  // legitimate (e.g., grm thresholds misordered)
      }
    }
  }

  UNPROTECT(1);
  return outsxp;
}

static SEXP
rpf_logprob_wrapper(SEXP r_spec, SEXP r_param, SEXP r_theta)
{
  if (length(r_spec) < RPF_ISpecCount)
    error("Item spec must be of length %d, not %d", RPF_ISpecCount, length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= librpf_numModels)
    error("Item model %d out of range", id);

  int numSpec = (*librpf_model[id].numSpec)(spec);
  if (length(r_spec) < numSpec)
    error("Item spec must be of length %d, not %d", numSpec, length(r_spec));
    
  int numParam = (*librpf_model[id].numParam)(spec);
  if (length(r_param) < numParam)
    error("Item has %d parameters, only %d given", numParam, length(r_param));

  const int numOutcomes = spec[RPF_ISpecOutcomes];
  const int dims = spec[RPF_ISpecDims];

  int numPeople;
  int numAbilities = 1;
  if (dims == 1) {
    numPeople = length(r_theta);
  } else {
    getMatrixDims(r_theta, &numAbilities, &numPeople);
  }

  SEXP outsxp;
  PROTECT(outsxp = allocMatrix(REALSXP, numOutcomes, numPeople));
  double *out = REAL(outsxp);
  double *theta = REAL(r_theta);
  double *param = REAL(r_param);
    
  for (int px=0; px < numPeople; px++) {
    double thBuf[dims];
    if (!unpack_theta(dims, param, numAbilities, theta + px*numAbilities, thBuf)) {
	for (int ox=0; ox < numOutcomes; ox++) {
	  out[px*numOutcomes + ox] = NA_REAL;
	}
	continue;
    }
    (*librpf_model[id].logprob)(spec, param, thBuf, out+px*numOutcomes);
    for (int ox=0; ox < numOutcomes; ox++) {
      double prob = out[px*numOutcomes + ox];
      if (!isfinite(prob)) {
	out[px*numOutcomes + ox] = NA_REAL;  // legitimate (e.g., grm thresholds misordered)
      }
    }
  }

  UNPROTECT(1);
  return outsxp;
}

static SEXP
rpf_dLL_wrapper(SEXP r_spec, SEXP r_param,
		  SEXP r_where, SEXP r_weight)
{
  if (length(r_spec) < RPF_ISpecCount)
    error("Item spec must be of length %d, not %d", RPF_ISpecCount, length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= librpf_numModels)
    error("Item model %d out of range", id);

  int numSpec = (*librpf_model[id].numSpec)(spec);
  if (length(r_spec) < numSpec)
    error("Item spec must be of length %d, not %d", numSpec, length(r_spec));
    
  int numParam = (*librpf_model[id].numParam)(spec);
  if (length(r_param) < numParam)
    error("Item has %d parameters, only %d given", numParam, length(r_param));

  int dims = spec[RPF_ISpecDims];
  if (length(r_where) != dims)
    error("Item has %d dimensions, but where is of length %d",
	  dims, length(r_where));

  int outcomes = spec[RPF_ISpecOutcomes];
  if (length(r_weight) != outcomes)
    error("Item has %d outcomes, but weight is of length %d",
	  outcomes, length(r_weight));

  const int numDeriv = numParam + numParam*(numParam+1)/2;
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, numDeriv));
  memset(REAL(ret), 0, sizeof(double) * numDeriv);
  (*librpf_model[id].dLL1)(spec, REAL(r_param),
			    REAL(r_where), REAL(r_weight), REAL(ret));
  for (int px=0; px < numDeriv; px++) {
    if (!isfinite(REAL(ret)[px])) error("Deriv %d not finite at step 1", px);
  }
  (*librpf_model[id].dLL2)(spec, REAL(r_param), REAL(ret));
  for (int px=0; px < numDeriv; px++) {
    if (!isfinite(REAL(ret)[px])) error("Deriv %d not finite at step 2", px);
  }
  UNPROTECT(1);
  return ret;
}

static SEXP
rpf_dTheta_wrapper(SEXP r_spec, SEXP r_param, SEXP r_where, SEXP r_dir)
{
  if (length(r_spec) < RPF_ISpecCount)
    error("Item spec must be of length %d, not %d", RPF_ISpecCount, length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= librpf_numModels)
    error("Item model %d out of range", id);

  int numSpec = (*librpf_model[id].numSpec)(spec);
  if (length(r_spec) < numSpec)
    error("Item spec must be of length %d, not %d", numSpec, length(r_spec));
    
  int numParam = (*librpf_model[id].numParam)(spec);
  if (length(r_param) < numParam)
    error("Item has %d parameters, only %d given", numParam, length(r_param));

  int dims = spec[RPF_ISpecDims];
  if (length(r_dir) != dims)
    error("Item has %d dimensions, but dir is of length %d",
	  dims, length(r_dir));
  if (length(r_where) != dims)
    error("Item has %d dimensions, but where is of length %d",
	  dims, length(r_where));

  SEXP ret, names;
  PROTECT(ret = allocVector(VECSXP, 2));
  PROTECT(names = allocVector(STRSXP, 2));

  int outcomes = spec[RPF_ISpecOutcomes];
  SEXP grad, hess;
  PROTECT(grad = allocVector(REALSXP, outcomes));
  PROTECT(hess = allocVector(REALSXP, outcomes));
  memset(REAL(grad), 0, sizeof(double) * outcomes);
  memset(REAL(hess), 0, sizeof(double) * outcomes);
  (*librpf_model[id].dTheta)(spec, REAL(r_param), REAL(r_where), REAL(r_dir),
			     REAL(grad), REAL(hess));
  SET_VECTOR_ELT(ret, 0, grad);
  SET_VECTOR_ELT(ret, 1, hess);
  SET_STRING_ELT(names, 0, mkChar("gradient"));
  SET_STRING_ELT(names, 1, mkChar("hessian"));
  namesgets(ret, names);
  UNPROTECT(4);
  return ret;
}

static SEXP
rpf_rescale_wrapper(SEXP r_spec, SEXP r_param, SEXP r_mean, SEXP r_cov)
{
  if (length(r_spec) < RPF_ISpecCount)
    error("Item spec must be of length %d, not %d", RPF_ISpecCount, length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= librpf_numModels)
    error("Item model %d out of range", id);

  int numSpec = (*librpf_model[id].numSpec)(spec);
  if (length(r_spec) < numSpec)
    error("Item spec must be of length %d, not %d", numSpec, length(r_spec));
    
  int numParam = (*librpf_model[id].numParam)(spec);
  if (length(r_param) < numParam)
    error("Item has %d parameters, only %d given", numParam, length(r_param));

  int dims = spec[RPF_ISpecDims];
  if (length(r_mean) != dims)
    error("Item has %d dimensions, but mean is of length %d",
	  dims, length(r_mean));

  int cov_rows, cov_cols;
  getMatrixDims(r_cov, &cov_rows, &cov_cols);
  if (cov_rows != dims || cov_rows != dims)
    error("Item has %d dimensions, but cov is %dx%d",
	  dims, cov_rows, cov_cols);

  int mask[numParam];
  memset(mask, 0, sizeof(*mask) * numParam);

  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, numParam));
  memcpy(REAL(ret), REAL(r_param), sizeof(double) * numParam);
  (*librpf_model[id].rescale)(spec, REAL(ret), mask,
			      REAL(r_mean), REAL(r_cov));
  UNPROTECT(1);
  return ret;
}

static R_CallMethodDef flist[] = {
  {"get_model_names", (DL_FUNC) get_model_names, 1},
  {"rpf_numSpec_wrapper", (DL_FUNC) rpf_numSpec_wrapper, 1},
  {"rpf_numParam_wrapper", (DL_FUNC) rpf_numParam_wrapper, 1},
  {"rpf_paramInfo_wrapper", (DL_FUNC) rpf_paramInfo_wrapper, 2},
  {"rpf_prob_wrapper", (DL_FUNC) rpf_prob_wrapper, 3},
  {"rpf_logprob_wrapper", (DL_FUNC) rpf_logprob_wrapper, 3},
  {"rpf_dLL_wrapper", (DL_FUNC) rpf_dLL_wrapper, 4},
  {"rpf_dTheta_wrapper", (DL_FUNC) rpf_dTheta_wrapper, 4},
  {"rpf_rescale_wrapper", (DL_FUNC) rpf_rescale_wrapper, 4},
  {"orlando_thissen_2000_wrapper", (DL_FUNC) orlando_thissen_2000, 5},
  {"kang_chen_2007_wrapper", (DL_FUNC) kang_chen_2007_wrapper, 2},
  {"sumscore_observed", (DL_FUNC) sumscore_observed, 4},
  {"ordinal_gamma_wrapper", (DL_FUNC) gamma_cor, 1},
  {NULL, NULL, 0}
};

static void
get_librpf_models(int *version, int *numModels, const struct rpf **model)
{
  *version = 8;
  *numModels = librpf_numModels;
  *model = librpf_model;
}

void R_init_rpf(DllInfo *info) {
  R_registerRoutines(info, NULL, flist, NULL, NULL);
  R_RegisterCCallable("rpf", "get_librpf_model_GPL", (DL_FUNC) get_librpf_models);
}

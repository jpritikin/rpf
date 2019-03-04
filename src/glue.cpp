#include "rpf.h"

SEXP MxRList::asR()
{
	// detect duplicate keys? TODO
	SEXP names, ans;
	int len = size();
	Rf_protect(names = Rf_allocVector(STRSXP, len));
	Rf_protect(ans = Rf_allocVector(VECSXP, len));
	for (int lx=0; lx < len; ++lx) {
		const char *p1 = (*this)[lx].first;
		SEXP p2 = (*this)[lx].second;
		if (!p1 || !p2) Rf_error("Attempt to return NULL pointer to R");
		SET_STRING_ELT(names, lx, Rf_mkChar(p1));
		SET_VECTOR_ELT(ans,   lx, p2);
	}
	Rf_namesgets(ans, names);
	return ans;
}

static SEXP
get_model_names(SEXP name)
{
  const char *target = CHAR(STRING_ELT(name, 0));
  SEXP ret;
  Rf_protect(ret = Rf_allocVector(REALSXP, 1));
  REAL(ret)[0] = NA_REAL;
  for (int sx=0; sx < Glibrpf_numModels; sx++) {
    if (strcmp(Glibrpf_model[sx].name, target) == 0) {
      REAL(ret)[0] = sx;
    }
  }
  UNPROTECT(1);
  return ret;
}

void getMatrixDims(SEXP r_theta, int *rows, int *cols)
{
    SEXP matrixDims;
    Rf_protect(matrixDims = Rf_getAttrib(r_theta, R_DimSymbol));
    int *dimList = INTEGER(matrixDims);
    *rows = dimList[0];
    *cols = dimList[1];
    UNPROTECT(1);
}

static SEXP
rpf_numSpec_wrapper(SEXP r_spec)
{
  if (Rf_length(r_spec) < RPF_ISpecCount)
    Rf_error("Item spec must be of length %d, not %d", RPF_ISpecCount, Rf_length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= Glibrpf_numModels)
    Rf_error("Item model %d out of range", id);

  int numSpec = (*Glibrpf_model[id].numSpec)(spec);

  SEXP ret;
  Rf_protect(ret = Rf_allocVector(INTSXP, 1));
  INTEGER(ret)[0] = numSpec;
  UNPROTECT(1);

  return ret;
}

static SEXP
rpf_numParam_wrapper(SEXP r_spec)
{
  if (Rf_length(r_spec) < RPF_ISpecCount)
    Rf_error("Item spec must be of length %d, not %d", RPF_ISpecCount, Rf_length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= Glibrpf_numModels)
    Rf_error("Item model %d out of range", id);

  int numParam = (*Glibrpf_model[id].numParam)(spec);

  SEXP ret;
  Rf_protect(ret = Rf_allocVector(INTSXP, 1));
  INTEGER(ret)[0] = numParam;
  UNPROTECT(1);

  return ret;
}

static SEXP
rpf_paramInfo_wrapper(SEXP r_spec, SEXP r_paramNum)
{
  if (Rf_length(r_spec) < RPF_ISpecCount)
    Rf_error("Item spec must be of length %d, not %d", RPF_ISpecCount, Rf_length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= Glibrpf_numModels)
    Rf_error("Item model %d out of range", id);

  int pnum = Rf_asInteger(r_paramNum);
  int numParam = (*Glibrpf_model[id].numParam)(spec);
  if (pnum < 0 || pnum >= numParam) Rf_error("Item model %d has %d parameters", id, numParam);

  const char *type;
  double upper, lower;
  (*Glibrpf_model[id].paramInfo)(spec, pnum, &type, &upper, &lower);

  int len = 3;
  SEXP names, ans;
  Rf_protect(names = Rf_allocVector(STRSXP, len));
  Rf_protect(ans = Rf_allocVector(VECSXP, len));
  int lx = 0;
  SET_STRING_ELT(names, lx, Rf_mkChar("type"));
  SET_VECTOR_ELT(ans,   lx, Rf_ScalarString(Rf_mkChar(type)));
  SET_STRING_ELT(names, ++lx, Rf_mkChar("upper"));
  SET_VECTOR_ELT(ans,   lx, Rf_ScalarReal(std::isfinite(upper)? upper : NA_REAL));
  SET_STRING_ELT(names, ++lx, Rf_mkChar("lower"));
  SET_VECTOR_ELT(ans,   lx, Rf_ScalarReal(std::isfinite(lower)? lower : NA_REAL));
  Rf_namesgets(ans, names);
  UNPROTECT(2);

  return ans;
}

int unpack_theta(int dims, double *param, int numAbilities, double *theta, double *out)
{
  if (numAbilities == dims) {
    for (int dx=0; dx < dims; ++dx) {
      double th = theta[dx];
      if (!std::isfinite(th)) return 0;
      out[dx] = th;
    }
  } else {
    int ax = 0;
    for (int dx=0; dx < dims; ++dx) {
      if (param[dx] == 0) continue;
      double th = theta[ax]; // could read uninitialized memory, but we detect below
      if (!std::isfinite(th)) return 0;
      out[dx] = th;
      ++ax;
    }
    if (ax != numAbilities) {
      Rf_error("Item has %d nonzero dims but given %d abilities", ax, numAbilities);
    }
  }
  return 1;
}

static SEXP
rpf_prob_wrapper(SEXP r_spec, SEXP r_param, SEXP r_theta)
{
  if (Rf_length(r_spec) < RPF_ISpecCount)
    Rf_error("Item spec must be of length %d, not %d", RPF_ISpecCount, Rf_length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= Glibrpf_numModels)
    Rf_error("Item model %d out of range", id);

  int numSpec = (*Glibrpf_model[id].numSpec)(spec);
  if (Rf_length(r_spec) < numSpec)
    Rf_error("Item spec must be of length %d, not %d", numSpec, Rf_length(r_spec));
    
  int numParam = (*Glibrpf_model[id].numParam)(spec);
  if (Rf_length(r_param) < numParam)
    Rf_error("Item has %d parameters, only %d given", numParam, Rf_length(r_param));

  const int numOutcomes = spec[RPF_ISpecOutcomes];
  const int dims = spec[RPF_ISpecDims];
  double *param = REAL(r_param);
  int numPeople = 1;
  int numAbilities = 1;
  if (dims == 0) {
	  if (Rf_length(r_theta)) numPeople = Rf_length(r_theta);
  } else if (dims == 1) {
    numPeople = Rf_length(r_theta);
  } else {
    getMatrixDims(r_theta, &numAbilities, &numPeople);
  }

  SEXP outsxp;
  Rf_protect(outsxp = Rf_allocMatrix(REALSXP, numOutcomes, numPeople));
  double *out = REAL(outsxp);
  double *theta = NULL;
  if (dims) theta = REAL(r_theta);
    
  Eigen::VectorXd thBuf(dims);
  for (int px=0; px < numPeople; px++) {
	  if (dims && !unpack_theta(dims, param, numAbilities, theta + px*numAbilities, thBuf.data())) {
		  for (int ox=0; ox < numOutcomes; ox++) {
			  out[px*numOutcomes + ox] = NA_REAL;
		  }
		  continue;
	  }
	  (*Glibrpf_model[id].prob)(spec, param, thBuf.data(), out+px*numOutcomes);
    for (int ox=0; ox < numOutcomes; ox++) {
      double prob = out[px*numOutcomes + ox];
      if (!std::isfinite(prob)) {
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
  if (Rf_length(r_spec) < RPF_ISpecCount)
    Rf_error("Item spec must be of length %d, not %d", RPF_ISpecCount, Rf_length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= Glibrpf_numModels)
    Rf_error("Item model %d out of range", id);

  int numSpec = (*Glibrpf_model[id].numSpec)(spec);
  if (Rf_length(r_spec) < numSpec)
    Rf_error("Item spec must be of length %d, not %d", numSpec, Rf_length(r_spec));
    
  int numParam = (*Glibrpf_model[id].numParam)(spec);
  if (Rf_length(r_param) < numParam)
    Rf_error("Item has %d parameters, only %d given", numParam, Rf_length(r_param));

  const int numOutcomes = spec[RPF_ISpecOutcomes];
  const int dims = spec[RPF_ISpecDims];
  int numPeople = 1;
  double *param = REAL(r_param);
  int numAbilities = 1;
  if (dims == 0) {
	  if (Rf_length(r_theta)) numPeople = Rf_length(r_theta);
  } else if (dims == 1) {
    numPeople = Rf_length(r_theta);
  } else {
    getMatrixDims(r_theta, &numAbilities, &numPeople);
  }

  SEXP outsxp;
  Rf_protect(outsxp = Rf_allocMatrix(REALSXP, numOutcomes, numPeople));
  double *out = REAL(outsxp);
  double *theta = NULL;
  if (dims) theta = REAL(r_theta);
    
  Eigen::VectorXd thBuf(dims);
  for (int px=0; px < numPeople; px++) {
	  if (dims && !unpack_theta(dims, param, numAbilities, theta + px*numAbilities, thBuf.data())) {
		  for (int ox=0; ox < numOutcomes; ox++) {
			  out[px*numOutcomes + ox] = NA_REAL;
		  }
		  continue;
	  }
	  (*Glibrpf_model[id].logprob)(spec, param, thBuf.data(), out+px*numOutcomes);
    for (int ox=0; ox < numOutcomes; ox++) {
      double prob = out[px*numOutcomes + ox];
      if (!std::isfinite(prob)) {
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
  if (Rf_length(r_spec) < RPF_ISpecCount)
    Rf_error("Item spec must be of length %d, not %d", RPF_ISpecCount, Rf_length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= Glibrpf_numModels)
    Rf_error("Item model %d out of range", id);

  int numSpec = (*Glibrpf_model[id].numSpec)(spec);
  if (Rf_length(r_spec) < numSpec)
    Rf_error("Item spec must be of length %d, not %d", numSpec, Rf_length(r_spec));
    
  int numParam = (*Glibrpf_model[id].numParam)(spec);
  if (Rf_length(r_param) < numParam)
    Rf_error("Item has %d parameters, only %d given", numParam, Rf_length(r_param));

  int dims = spec[RPF_ISpecDims];
  if (Rf_length(r_where) != dims)
    Rf_error("Item has %d dimensions, but where is of length %d",
	  dims, Rf_length(r_where));

  int outcomes = spec[RPF_ISpecOutcomes];
  if (Rf_length(r_weight) != outcomes)
    Rf_error("Item has %d outcomes, but weight is of length %d",
	  outcomes, Rf_length(r_weight));

  double *where = NULL;
  if (dims) where = REAL(r_where);

  const int numDeriv = numParam + numParam*(numParam+1)/2;
  SEXP ret;
  Rf_protect(ret = Rf_allocVector(REALSXP, numDeriv));
  memset(REAL(ret), 0, sizeof(double) * numDeriv);
  (*Glibrpf_model[id].dLL1)(spec, REAL(r_param),
			    where, REAL(r_weight), REAL(ret));
  for (int px=0; px < numDeriv; px++) {
    if (!std::isfinite(REAL(ret)[px])) Rf_error("Deriv %d not finite at step 1", px);
  }
  (*Glibrpf_model[id].dLL2)(spec, REAL(r_param), REAL(ret));
  for (int px=0; px < numDeriv; px++) {
	  //if (!std::isfinite(REAL(ret)[px])) Rf_error("Deriv %d not finite at step 2", px);
  }
  UNPROTECT(1);
  return ret;
}

static SEXP
rpf_dTheta_wrapper(SEXP r_spec, SEXP r_param, SEXP r_where, SEXP r_dir)
{
  if (Rf_length(r_spec) < RPF_ISpecCount)
    Rf_error("Item spec must be of length %d, not %d", RPF_ISpecCount, Rf_length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= Glibrpf_numModels)
    Rf_error("Item model %d out of range", id);

  int numSpec = (*Glibrpf_model[id].numSpec)(spec);
  if (Rf_length(r_spec) < numSpec)
    Rf_error("Item spec must be of length %d, not %d", numSpec, Rf_length(r_spec));
    
  int numParam = (*Glibrpf_model[id].numParam)(spec);
  if (Rf_length(r_param) < numParam)
    Rf_error("Item has %d parameters, only %d given", numParam, Rf_length(r_param));

  int dims = spec[RPF_ISpecDims];
  if (dims == 0) Rf_error("Item has no factors");
  if (Rf_length(r_dir) != dims)
    Rf_error("Item has %d dimensions, but dir is of length %d",
	  dims, Rf_length(r_dir));
  if (Rf_length(r_where) != dims)
    Rf_error("Item has %d dimensions, but where is of length %d",
	  dims, Rf_length(r_where));

  SEXP ret, names;
  Rf_protect(ret = Rf_allocVector(VECSXP, 2));
  Rf_protect(names = Rf_allocVector(STRSXP, 2));

  int outcomes = spec[RPF_ISpecOutcomes];
  SEXP grad, hess;
  Rf_protect(grad = Rf_allocVector(REALSXP, outcomes));
  Rf_protect(hess = Rf_allocVector(REALSXP, outcomes));
  memset(REAL(grad), 0, sizeof(double) * outcomes);
  memset(REAL(hess), 0, sizeof(double) * outcomes);
  (*Glibrpf_model[id].dTheta)(spec, REAL(r_param), REAL(r_where), REAL(r_dir),
			     REAL(grad), REAL(hess));
  SET_VECTOR_ELT(ret, 0, grad);
  SET_VECTOR_ELT(ret, 1, hess);
  SET_STRING_ELT(names, 0, Rf_mkChar("gradient"));
  SET_STRING_ELT(names, 1, Rf_mkChar("hessian"));
  Rf_namesgets(ret, names);
  UNPROTECT(4);
  return ret;
}

static SEXP
rpf_rescale_wrapper(SEXP r_spec, SEXP r_param, SEXP r_mean, SEXP r_cov)
{
  if (Rf_length(r_spec) < RPF_ISpecCount)
    Rf_error("Item spec must be of length %d, not %d", RPF_ISpecCount, Rf_length(r_spec));

  double *spec = REAL(r_spec);

  int id = spec[RPF_ISpecID];
  if (id < 0 || id >= Glibrpf_numModels)
    Rf_error("Item model %d out of range", id);

  int numSpec = (*Glibrpf_model[id].numSpec)(spec);
  if (Rf_length(r_spec) < numSpec)
    Rf_error("Item spec must be of length %d, not %d", numSpec, Rf_length(r_spec));
    
  int numParam = (*Glibrpf_model[id].numParam)(spec);
  if (Rf_length(r_param) < numParam)
    Rf_error("Item has %d parameters, only %d given", numParam, Rf_length(r_param));

  int dims = spec[RPF_ISpecDims];
  if (dims == 0) Rf_error("Item has no factors");
  if (Rf_length(r_mean) != dims)
    Rf_error("Item has %d dimensions, but mean is of length %d",
	  dims, Rf_length(r_mean));

  int cov_rows, cov_cols;
  getMatrixDims(r_cov, &cov_rows, &cov_cols);
  if (cov_rows != dims || cov_rows != dims)
    Rf_error("Item has %d dimensions, but cov is %dx%d",
	  dims, cov_rows, cov_cols);

  Eigen::VectorXi mask(numParam);
  mask.setZero();

  SEXP ret;
  Rf_protect(ret = Rf_allocVector(REALSXP, numParam));
  memcpy(REAL(ret), REAL(r_param), sizeof(double) * numParam);
  (*Glibrpf_model[id].rescale)(spec, REAL(ret), mask.data(),
			      REAL(r_mean), REAL(r_cov));
  UNPROTECT(1);
  return ret;
}

static SEXP has_openmp()
{
#if defined(_OPENMP)
	bool opm = true;
#else
	bool opm = false;
#endif
	return Rf_ScalarLogical(opm);
}

int GlobalNumberOfCores = 1;

static SEXP setNumberOfCores(SEXP num)
{
#if defined(_OPENMP)
	GlobalNumberOfCores = Rf_asInteger(num);
#endif
	return num;
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
  {"collapse_wrapper", (DL_FUNC) collapse_wrapper, 2},
  {"ordinal_gamma_wrapper", (DL_FUNC) gamma_cor, 1},
  {"ssEAP_wrapper", (DL_FUNC) sumscoreEAP, 6},
  {"ot2000_wrapper", (DL_FUNC) ot2000_wrapper, 7},
  {"crosstabTest_wrapper", (DL_FUNC) crosstabTest, 3},
  {"pairwiseExpected_wrapper", (DL_FUNC) pairwiseExpected, 5},
  {"observedSumScore_wrapper", (DL_FUNC) observedSumScore, 2},
  {"itemOutcomeBySumScore_wrapper", (DL_FUNC) itemOutcomeBySumScore, 3},
  {"findIdenticalRowsData", (DL_FUNC) findIdenticalRowsData, 5},
  {"CaiHansen2012_wrapper", (DL_FUNC) CaiHansen2012, 3},
  {"eap_wrapper", (DL_FUNC) eap_wrapper, 1},
  {"hasOpenMP_wrapper", (DL_FUNC) has_openmp, 0},
  {"setNumberOfCores", (DL_FUNC) setNumberOfCores, 1},
  {"fast_tableWithWeights", (DL_FUNC) fast_tableWithWeights, 3},
  {NULL, NULL, 0}
};

extern const struct rpf librpf_model[];
extern const int librpf_numModels;

static void
get_librpf_models(int version, int *numModels, const struct rpf **model)
{
	if (version != LIBIFA_RPF_API_VERSION) Rf_error("LIBIFA_RPF binary API version mismatch");
  *numModels = librpf_numModels;
  *model = librpf_model;
}

const struct rpf *Glibrpf_model;
int Glibrpf_numModels;

extern "C" void R_init_rpf(DllInfo *info) {
	R_registerRoutines(info, NULL, flist, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
	R_RegisterCCallable("rpf", "get_librpf_model_GPL", (DL_FUNC) get_librpf_models);
	Glibrpf_numModels = librpf_numModels;
	Glibrpf_model = librpf_model;
}

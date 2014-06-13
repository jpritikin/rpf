// All source files should include this single header file. The idea here
// is to make it easy to switch to using precompiled headers. Currently,
// the full compile is fast enough that I haven't bothered.

#ifdef EIGEN_WORLD_VERSION
#error "rpf.h must be included before Eigen"
#endif

#define EIGEN_NO_DEBUG 1
#define EIGEN_DONT_PARALLELIZE

#ifdef DEBUG
#define EIGEN_INITIALIZE_MATRICES_BY_NAN
#undef EIGEN_NO_DEBUG
#endif

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <vector>
#include "Eigen/Core"
#include "../inst/include/libifa-rpf.h"
#include "dmvnorm.h"
#include "ba81quad.h"

class omxManageProtectInsanity {
	PROTECT_INDEX initialpix;
 public:
	omxManageProtectInsanity() {
		PROTECT_WITH_INDEX(R_NilValue, &initialpix);
		UNPROTECT(1);
	}
	PROTECT_INDEX getDepth() {
		PROTECT_INDEX pix;
		PROTECT_WITH_INDEX(R_NilValue, &pix);
		PROTECT_INDEX diff = pix - initialpix;
		UNPROTECT(1);
		return diff;
	}
	~omxManageProtectInsanity() {
		UNPROTECT(getDepth());
	}
};

typedef std::vector< std::pair<const char *, SEXP> > MxRListBase;
class MxRList : private MxRListBase {
 public:
	size_t size() const { return MxRListBase::size(); }
	SEXP asR();
	void add(const char *key, SEXP val) {
		Rf_protect(val);
		push_back(std::make_pair(key, val));
	};
};

static inline bool strEQ(const char *s1, const char *s2) { return strcmp(s1,s2)==0; }

void getMatrixDims(SEXP r_theta, int *rows, int *cols);
SEXP orlando_thissen_2000(SEXP r_spec, SEXP r_param, SEXP r_item, SEXP r_observed, SEXP r_quad);
SEXP collapse_wrapper(SEXP r_observed_orig, SEXP r_expected_orig);
SEXP gamma_cor(SEXP r_mat);
SEXP sumscoreEAP(SEXP robj, SEXP Rwidth, SEXP Rpts, SEXP Rmask);
SEXP ot2000_wrapper(SEXP robj, SEXP Ritem, SEXP Rwidth, SEXP Rpts, SEXP Ralter, SEXP Rmask);
SEXP crosstabTest(SEXP Robserved, SEXP Rexpected, SEXP Rtrials);
SEXP pairwiseExpected(SEXP robj, SEXP Rwidth, SEXP Rpts, SEXP Ritems);
SEXP observedSumScore(SEXP Rgrp, SEXP Rmask);
SEXP itemOutcomeBySumScore(SEXP Rgrp, SEXP Rmask, SEXP Rinterest);

static inline int triangleLoc1(int diag)
{
	return (diag) * (diag+1) / 2;   // 0 1 3 6 10 15 ..
}

static inline int triangleLoc0(int diag)
{
	return triangleLoc1(diag+1) - 1;  // 0 2 5 9 14 ..
}

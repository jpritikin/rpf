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

#include <Rcpp.h>
using namespace Rcpp;

#include "Eigen/Core"

#if defined(_OPENMP)
#include <omp.h>
#else
static inline int omp_get_thread_num() { return 0; }
#endif

static inline int triangleLoc1(int diag)
{
	return (diag) * (diag+1) / 2;   // 0 1 3 6 10 15 ..
}

static inline int triangleLoc0(int diag)
{
	return triangleLoc1(diag+1) - 1;  // 0 2 5 9 14 ..
}

#include "../inst/include/libifa-rpf.h"
#include "dmvnorm.h"
#include "ba81quad.h"

extern int GlobalNumberOfCores;

class ProtectAutoBalanceDoodad {
	PROTECT_INDEX initialpix;
 public:
	ProtectAutoBalanceDoodad() {
		R_ProtectWithIndex(R_NilValue, &initialpix);
		Rf_unprotect(1);
	}
	PROTECT_INDEX getDepth() {
		PROTECT_INDEX pix;
		R_ProtectWithIndex(R_NilValue, &pix);
		PROTECT_INDEX diff = pix - initialpix;
		Rf_unprotect(1);
		return diff;
	}
	~ProtectAutoBalanceDoodad() {
		Rf_unprotect(getDepth());
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

static inline void
pda(const double *ar, int rows, int cols) {   // column major order
	for (int rx=0; rx < rows; rx++) {
		for (int cx=0; cx < cols; cx++) {
			Rprintf("%.6g, ", ar[cx * rows + rx]);
		}
		Rprintf("\n");
	}
}

static inline void string_to_try_Rf_error( const std::string& str )
{
	Rf_error("%s", str.c_str());
}

static inline void exception_to_try_Rf_error( const std::exception& ex )
{
	string_to_try_Rf_error(ex.what());
}

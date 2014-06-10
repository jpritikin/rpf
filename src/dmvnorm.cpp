#include "rpf.h"

struct Matrix {
	int rows;
	int cols;
	double *t;

	Matrix() {}
	Matrix(double *_t, int _r, int _c) : rows(_r), cols(_c), t(_t) {}
};

int MatrixSolve(Matrix mat1, Matrix mat2, bool identity)
{
	if (mat1.rows != mat1.cols ||
	    mat2.rows != mat2.cols ||
	    mat1.rows != mat2.rows) Rf_error("Not conformable");
	const int dim = mat1.rows;

	Eigen::VectorXd pad(dim * dim);
	memcpy(pad.data(), mat1.t, sizeof(double) * dim * dim);

	if (identity) {
		for (int rx=0; rx < dim; rx++) {
			for (int cx=0; cx < dim; cx++) {
				mat2.t[rx * dim + cx] = rx==cx? 1 : 0;
			}
		}
	}
  
	std::vector<int> ipiv(dim);
	int info;
	F77_CALL(dgesv)(&dim, &dim, pad.data(), &dim, ipiv.data(), mat2.t, &dim, &info);
	if (info < 0) {
		Rf_error("Arg %d is invalid", -info);
	}
	return info;
}

static const int ERROR_LEN = 80;

static double
_mahalanobis(char *err, int dim, double *loc, double *center, double *origCov)
{
	std::vector<double> cloc(dim);
	for (int dx=0; dx < dim; dx++) {
		cloc[dx] = loc[dx] - center[dx];
	}

	Matrix covMat(origCov, dim, dim);
	std::vector<double> icov(dim * dim);
	Matrix icovMat(icov.data(), dim, dim);
	int info = MatrixSolve(covMat, icovMat, true); // can optimize for symmetry TODO
	if (info) {
		snprintf(err, ERROR_LEN, "Sigma is singular and cannot be inverted");
		return nan("Rf_error");
	}

	std::vector<double> half(dim);
	char trans='n';
	double alpha=1;
	double beta=0;
	int inc=1;
	F77_CALL(dgemv)(&trans, &dim, &dim, &alpha, icov.data(), &dim, cloc.data(), &inc, &beta, half.data(), &inc);

	double got=0;
	for (int dx=0; dx < dim; dx++) got += half[dx] * cloc[dx];
	return got;
}

static double
mahalanobis(int dim, double *loc, double *center, double *origCov)
{
	char err[ERROR_LEN];
	err[0] = 0;
	double ret = _mahalanobis(err, dim, loc, center, origCov);
	if (err[0]) Rf_error("%s", err);
	return ret;
}

static double
_dmvnorm(char *err, int dim, double *loc, double *mean, double *origSigma)
{
	double dist = mahalanobis(dim, loc, mean, origSigma);

	std::vector<double> sigma(dim * dim);
	memcpy(sigma.data(), origSigma, sizeof(double) * dim * dim);

	char jobz = 'N';
	char range = 'A';
	char uplo = 'U';
	double vunused;
	int iunused;
	double abstol = 0;
	int m;
	Eigen::VectorXd w(dim);
	Eigen::VectorXd Z(dim);
	int ldz=1;
	Eigen::VectorXi isuppz(2*dim);
	int lwork = -1;
	double optlWork;
	int optliWork;
	int liwork = -1;
	int info;

	F77_CALL(dsyevr)(&jobz, &range, &uplo,
			 &dim, sigma.data(), &dim,
			 &vunused, &vunused,
			 &iunused, &iunused,
			 &abstol, &m, w.data(),
			 Z.data(), &ldz, isuppz.data(),
			 &optlWork, &lwork,
			 &optliWork, &liwork, &info);
	if (info != 0) {
		snprintf(err, ERROR_LEN, "dsyevr failed when requesting work space size");
		return nan("Rf_error");
	}

	lwork = optlWork;
	std::vector<double> work(lwork);
	liwork = optliWork;
	std::vector<int> iwork(liwork);

	F77_CALL(dsyevr)(&jobz, &range, &uplo, &dim, sigma.data(), &dim,
			 &vunused, &vunused, &iunused, &iunused, &abstol, &m, w.data(), Z.data(), &ldz, isuppz.data(),
			 work.data(), &lwork, iwork.data(), &liwork, &info);
	if (info < 0) {
		snprintf(err, ERROR_LEN, "Arg %d is invalid", -info);
		return nan("Rf_error");
	}
	if (info > 0) {
		snprintf(err, ERROR_LEN, "dsyevr: internal Rf_error");
		return nan("Rf_error");
	}
	if (m < dim) {
		snprintf(err, ERROR_LEN, "Sigma not of full rank");
		return nan("Rf_error");
	}

	for (int dx=0; dx < dim; dx++) dist += log(w[dx]);
	double got = -(dim * M_LN_SQRT_2PI*2 + dist)/2;
	return got;
}

double
dmvnorm(int dim, double *loc, double *mean, double *sigma)
{
	char err[ERROR_LEN];
	err[0] = 0;
	double ret = _dmvnorm(err, dim, loc, mean, sigma);
	if (err[0]) Rf_error("%s", err);
	return ret;
}

SEXP dmvnorm_wrapper(SEXP Rloc, SEXP Rmean, SEXP Rsigma)
{
	SEXP ret;
	Rf_protect(ret = Rf_allocVector(REALSXP, 1));
	REAL(ret)[0] = dmvnorm(Rf_length(Rloc), REAL(Rloc), REAL(Rmean), REAL(Rsigma));
	Rf_unprotect(1);
	return ret;
}

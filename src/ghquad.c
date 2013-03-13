/**
 * Copyright 2011 Alexander W Blocker <ablocker@gmail.com>
 * 
 * Copied from version 0.1-1 of the fastGHQuad on CRAN.
 * MIT license.
 *
 * Originally C++, modified to compile in C.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

static void buildHermiteJacobi( int n, double* D, double* E ) {
//
// Construct symmetric tridiagonal matrix similar to Jacobi matrix
// for Hermite polynomials
//
// On exit, D contains diagonal elements of said matrix;
// E contains subdiagonal elements.
//
// Need D of size n, E of size n-1
//
// Building matrix based on recursion relation for monic versions of Hermite
// polynomials:
//      p_n(x) = H_n(x) / 2^n
//      p_n+1(x) + (B_n-x)*p_n(x) + A_n*p_n-1(x) = 0
//      B_n = 0
//      A_n = n/2
// 
// Matrix similar to Jacobi (J) defined by:
//      J_i,i = B_i-1, i = 1, ..., n
//      J_i,i-1 = J_i-1,i = sqrt(A_i-1), i = 2, ..., n
//

    // Build diagonal
    int i;
    for (i=0; i<n; i++) {
        D[i]    = 0;
    }

    // Build sub/super-diagonal
    for (i=0; i<n-1; i++) {
        E[i]    = sqrt((i+1.0)/2);
    }
    
}

static void quadInfoGolubWelsch( int n, double* D, double* E, double mu0,
                          double* x, double* w ) {
//
// Compute weights & nodes for Gaussian quadrature using Golub-Welsch algorithm.
//
// First need to build symmetric tridiagonal matrix J similar to Jacobi for
// desired orthogonal polynomial (based on recurrence relation).
// 
// D contains the diagonal of this matrix J, and E contains the
// sub/super-diagonal.
//
// This routine finds the eigenvectors & values of the given J matrix.
//
// The eigenvalues correspond to the nodes for the desired quadrature rule
// (roots of the orthogonal polynomial).
//
// The eigenvectors can be used to compute the weights for the quadrature rule
// via:
//
//      w_j = mu0 * (v_j,1)^2
// 
// where mu0 = \int_a^b w(x) dx
// (integral over range of integration of weight function)
//
// and 
//
// v_j,1 is the first entry of the jth normalized (to unit length) eigenvector.
//
// On exit, x (length n) contains nodes for quadrature rule, and w (length n)
// contains weights for quadrature rule.
//
// Note that contents of D & E are destroyed on exit
//

    // Setup for eigenvalue computations
    char JOBZ   = 'V'; // Compute eigenvalues & vectors
    int INFO;

    // Initialize array for workspace
    double * WORK   = Realloc(NULL, 2*n-2, double);

    // Initialize array for eigenvectors
    double * Z      = Realloc(NULL, n*n, double);

    // Run eigen decomposition
    F77_NAME(dstev)(&JOBZ, &n, D, E,     // Job flag & input matrix
            Z, &n,              // Output array for eigenvectors & dim
            WORK, &INFO         // Workspace & info flag
            );

    // Setup x & w
    int i;
    for (i=0; i<n; i++) {
        x[i] = D[i];
        w[i] = mu0*Z[i*n]*Z[i*n];
    }

    // Deallocate temporary arrays
    Free(WORK);
    Free(Z);
}

static void gaussHermiteData( int n, double* x, double* w ) {
//
// Calculates nodes & weights for Gauss-Hermite integration of order n
//
// Need x & w of size n
//
// Using standard formulation (no generalizations or polynomial adjustment)
//
// Evaluations use Golub-Welsch algorithm; numerically stable for n>=100
//
    // Build Jacobi-similar symmetric tridiagonal matrix via diagonal &
    // sub-diagonal
    double * D, * E;
    D   = Realloc(NULL, n, double);
    E   = Realloc(NULL, n, double);
    //
    buildHermiteJacobi( n, D, E );

    // Get nodes & weights
    double mu0  = M_SQRT_PI;
    quadInfoGolubWelsch( n, D, E, mu0,
                          x, w );

    // Deallocate temporary objects
    Free(D);
    Free(E);
}

SEXP omxGaussHermiteData(SEXP points)
{
	int pts = asInteger(points);
	SEXP ans;
	SEXP names;
	SEXP where, area;

	PROTECT(names = allocVector(STRSXP, 2));
	SET_STRING_ELT(names, 0, mkChar("x"));
	SET_STRING_ELT(names, 1, mkChar("w"));

	PROTECT(ans = allocVector(VECSXP, 2));
	PROTECT(where = allocVector(REALSXP, pts));
	SET_VECTOR_ELT(ans, 0, where);
	PROTECT(area = allocVector(REALSXP, pts));
	SET_VECTOR_ELT(ans, 1, area);

	gaussHermiteData(pts, REAL(where), REAL(area));

	namesgets(ans, names);
	UNPROTECT(4);

	return ans;
}

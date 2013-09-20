#include "rpf.h"

static void
decodeLocation(long qx, const int dims, const int gridsize,
	       const double *restrict pts, double *restrict out)
{
  for (int dx=0; dx < dims; dx++) {
    out[dx] = pts[qx % gridsize];
    qx = qx / gridsize;
  }
}

static void
pda(const double *ar, int rows, int cols) {
	for (int rx=0; rx < rows; rx++) {
		for (int cx=0; cx < cols; cx++) {
			Rprintf("%.6g ", ar[cx * rows + rx]);
		}
		Rprintf("\n");
	}

}

static const double KANG_CHEN_MIN_EXPECTED = 1.0;  // customizable parameter?

static void find_smallcol(int rows, int cols, double *expected, int rx,
			  int *goodcol, int *smallcol)
{
  *goodcol = 0;
  *smallcol = -1;
  double smallest = -1;
  for (int cx=0; cx < cols; cx++) {
    double cell = expected[cx * rows + rx];
    if (cell != 0 && cell <= KANG_CHEN_MIN_EXPECTED) {
      if (smallest == -1 || smallest > cell) {
	smallest = cell;
	*smallcol = cx;
      }
    } else {
      *goodcol += 1;
    }
  }
}

static int kang_chen_2007_collapse1(int rows, int cols, int *observed, double *expected, int rx)
{
  int collapsed = 0;
  int smallcol;
  int goodcol;
  find_smallcol(rows, cols, expected, rx, &goodcol, &smallcol);
  if (smallcol == -1) return collapsed;

  if (goodcol > cols/2) {
    // merge left or right
    while (1) {
      double biggest = 0;
      int bigcol = -1;
      if (smallcol - 1 >= 0) {
	double cell = expected[(smallcol-1)*rows + rx];
	if (biggest < cell) {
	  biggest = cell;
	  bigcol = smallcol-1;
	}
      }
      if (smallcol + 1 < cols) {
	double cell = expected[(smallcol+1)*rows + rx];
	if (biggest < cell) {
	  biggest = cell;
	  bigcol = smallcol+1;
	}
      }
      if (bigcol==-1) error("Confused");
      //Rprintf("collapse col %d to %d on row %d\n", smallcol, bigcol, rx);
      expected[bigcol*rows + rx] += expected[smallcol*rows + rx];
      observed[bigcol*rows + rx] += observed[smallcol*rows + rx];
      expected[smallcol*rows + rx] = 0;
      observed[smallcol*rows + rx] = 0;
      ++collapsed;

      find_smallcol(rows, cols, expected, rx, &goodcol, &smallcol);
      if (smallcol == -1) break;
    }
  } else {
    // merge up or down
    int bigrow;
    if (rx < rows/2) {
      bigrow = rx+1;
    } else {
      bigrow = rx-1;
    }
    for (int cx=0; cx < cols; cx++) {
      double e1 = expected[cx*rows + rx];
      if (e1 == 0 || e1 >= KANG_CHEN_MIN_EXPECTED)
	continue;
      //Rprintf("collapse row %d to %d at col %d\n", rx, bigrow, cx);
      expected[cx*rows + bigrow] += expected[cx*rows + rx];
      observed[cx*rows + bigrow] += observed[cx*rows + rx];
      expected[cx*rows + rx] = 0;
      observed[cx*rows + rx] = 0;
      ++collapsed;
    }
  }
  return collapsed;
}

/*
 * This assumes that the expected counts are greater than 1 in the
 * middle of the table. This may be suboptimal for large tables with
 * isolated regions of low expected counts (is it possible?). It may
 * be necessary to search the whole table at every step for the lowest
 * expected count and fix them in that order instead of assuming
 * anything about where the low expected counts will appear.
 */
static int kang_chen_2007_collapse(int rows, int cols, int *observed, double *expected)
{
  //Rprintf("kang chen 2007\n");
  //pda(expected, rows, cols);
  int collapsed = 0;
  for (int rx=0; rx <= rows/2; rx++) {
    collapsed += kang_chen_2007_collapse1(rows, cols, observed, expected, rx);
    if (rx != ((rows-1) -rx)) {
      collapsed += kang_chen_2007_collapse1(rows, cols, observed, expected, (rows-1) -rx);
    }
  }
  for (int rx=0; rx < rows; rx++) {
    for (int cx=0; cx < cols; cx++) {
      double cell = expected[cx * rows + rx];
      if (cell != 0 && cell <= KANG_CHEN_MIN_EXPECTED) {
	pda(expected, rows, cols);
	error("Failed to collapse cells");
      }
    }
  }
  return collapsed;
}

SEXP kang_chen_2007_wrapper(SEXP r_observed_orig, SEXP r_expected_orig)
{
  int rows, cols;
  getMatrixDims(r_expected_orig, &rows, &cols);

  {
    int orows, ocols;
    getMatrixDims(r_observed_orig, &orows, &ocols);
    if (rows != orows || cols != ocols)
      error("Observed and expected matrices must have same dimensions");
  }

  SEXP r_observed, r_expected;
  PROTECT(r_observed = allocMatrix(INTSXP, rows, cols));
  PROTECT(r_expected = allocMatrix(REALSXP, rows, cols));

  int *observed = INTEGER(r_observed);
  double *expected = REAL(r_expected);
  memcpy(observed, INTEGER(r_observed_orig), sizeof(int) * rows * cols);
  memcpy(expected, REAL(r_expected_orig), sizeof(double) * rows * cols);

  int collapsed = kang_chen_2007_collapse(rows, cols, observed, expected);

  const int returnCount = 3;
  SEXP names, ans;
  PROTECT(names = allocVector(STRSXP, returnCount));
  PROTECT(ans = allocVector(VECSXP, returnCount));

  int ansC = -1;
  SET_STRING_ELT(names, ++ansC, mkChar("observed"));
  SET_VECTOR_ELT(ans,   ansC, r_observed);
  SET_STRING_ELT(names, ++ansC, mkChar("expected"));
  SET_VECTOR_ELT(ans,   ansC, r_expected);
  SET_STRING_ELT(names, ++ansC, mkChar("collapsed"));
  SET_VECTOR_ELT(ans,   ansC, ScalarInteger(collapsed));
  if (ansC != returnCount-1) error("Memory corruption");

  namesgets(ans, names);
  UNPROTECT(4);

  return ans;
}

static const int CRAZY_VERSION = 0;  // remove TODO

SEXP orlando_thissen_2000(SEXP r_spec, SEXP r_param, SEXP r_item, SEXP r_observed_orig, SEXP r_quad)
{
  int numSpec = length(r_spec);
  if (numSpec < 2) error("At least 2 items are needed");
  int maxParam;
  int numItems;
  getMatrixDims(r_param, &maxParam, &numItems);
  if (numItems != numSpec) {
    error("Number of items %d and item specifications %d need to match", numItems, numSpec);
  }

  double *spec[numItems];
  int contextMap[numItems];
  int maxDims = 1;
  int maxCorrect = 0;
  int nextContext = 0;
  int interest = asInteger(r_item) - 1;
  if (interest < 0 || interest >= numSpec) {
    error("Item of interest %d must be between 1 and %d", interest+1, numSpec);
  }
  for (int sx=0; sx < numSpec; sx++) {
    spec[sx] = REAL(VECTOR_ELT(r_spec, sx));
    int dims = spec[sx][RPF_ISpecDims];
    if (maxDims < dims)
      maxDims = dims;
    maxCorrect += spec[sx][RPF_ISpecOutcomes] - 1;

    if (sx != interest) contextMap[nextContext++] = sx;
  }
  contextMap[nextContext] = interest; // remove TODO

  double *param = REAL(r_param);
  int quad_size = length(VECTOR_ELT(r_quad, 0));
  double *quad_pts = REAL(VECTOR_ELT(r_quad, 0));
  double *quad_area = REAL(VECTOR_ELT(r_quad, 1));
  int totalQuadOrdinate = 1;
  for (int dx=0; dx < maxDims; dx++) { totalQuadOrdinate *= quad_size; }
  double *lk1 = Realloc(NULL, (1+maxCorrect) * totalQuadOrdinate, double);
  memset(lk1, 0, sizeof(double) * (1+maxCorrect) * totalQuadOrdinate); // TODO remove
  double *lk2 = Realloc(NULL, (1+maxCorrect) * totalQuadOrdinate, double);

  double *curSpec = spec[contextMap[0]];
  int curOutcomes = curSpec[RPF_ISpecOutcomes];
  int curMaxCorrect = curOutcomes - 1;
  int id = curSpec[RPF_ISpecID];
  for (int qx=0; qx < totalQuadOrdinate; qx++) {
    double where[maxDims];
    decodeLocation(qx, maxDims, quad_size, quad_pts, where);
    double oprob[curOutcomes];
    librpf_model[id].prob(curSpec, param + maxParam * contextMap[0], where, oprob);
    for (int ox=0; ox < curOutcomes; ox++) {
      double *lk = lk1 + ox * totalQuadOrdinate;
      lk[qx] = oprob[ox];
    }
  }
  //pda(lk1, totalQuadOrdinate, 1+maxCorrect);

  for (int ix=1; ix < numItems; ix++) {  // TODO don't need last layer
    memset(lk2, 0, sizeof(double) * (1+maxCorrect) * totalQuadOrdinate);
    curSpec = spec[contextMap[ix]];
    curOutcomes = curSpec[RPF_ISpecOutcomes];
    id = curSpec[RPF_ISpecID];
    for (int qx=0; qx < totalQuadOrdinate; qx++) {
      double where[maxDims];
      decodeLocation(qx, maxDims, quad_size, quad_pts, where);
      double oprob[curOutcomes];
      librpf_model[id].prob(curSpec, param + maxParam * contextMap[ix], where, oprob);
      for (int cx=0; cx <= curMaxCorrect; cx++) {
	double *oldlk = lk1 + cx * totalQuadOrdinate;
	for (int ox=0; ox < curOutcomes; ox++) {
	  double *lk = lk2 + (cx + ox) * totalQuadOrdinate;
	  lk[qx] += oldlk[qx] * oprob[ox];
	}
      }
    }
    { double *tmp = lk1; lk1 = lk2; lk2 = tmp; }
    curMaxCorrect += curOutcomes - 1;
    //pda(lk1, totalQuadOrdinate, 1+maxCorrect);
  }
  curMaxCorrect -= curOutcomes - 1;

  // lk1 = with item of interest
  // lk2 = without item of interest

  double denom[maxCorrect+1];
  memset(denom, 0, sizeof(*denom) * (maxCorrect+1));
  for (int cx=0; cx <= maxCorrect; cx++) {
    double *lk;
    if (CRAZY_VERSION) {
      lk = lk1 + cx * totalQuadOrdinate;
    } else {
      lk = lk2 + cx * totalQuadOrdinate;
    }
    for (int qx=0; qx < totalQuadOrdinate; qx++) {
      double area[maxDims];
      decodeLocation(qx, maxDims, quad_size, quad_area, area);
      double piece = lk[qx];
      for (int dx=0; dx < maxDims; dx++) piece *= area[dx];
      denom[cx] += piece;
    }
  }
  //pda(denom, curMaxCorrect+1, 1);

  int outRows;
  if (CRAZY_VERSION) {
    outRows = maxCorrect+1;
  } else {
    outRows = curMaxCorrect+1;
  }

  SEXP r_expected_orig;
  PROTECT(r_expected_orig = allocMatrix(REALSXP, outRows, curOutcomes));
  double *expected_orig = REAL(r_expected_orig);
  memset(expected_orig, 0, sizeof(double) * outRows * curOutcomes);

  curSpec = spec[interest];
  curOutcomes = curSpec[RPF_ISpecOutcomes];
  id = curSpec[RPF_ISpecID];

  for (int qx=0; qx < totalQuadOrdinate; qx++) {
    double where[maxDims];
    decodeLocation(qx, maxDims, quad_size, quad_pts, where);
    double area[maxDims];
    decodeLocation(qx, maxDims, quad_size, quad_area, area);

    for (int cx=0; cx < outRows; cx++) {
      double oprob[curOutcomes];
      librpf_model[id].prob(curSpec, param + maxParam * interest, where, oprob);

      for (int ox=0; ox < curOutcomes; ox++) {
	double *lk;
	if (CRAZY_VERSION) {
	  lk = lk2 + (cx-ox) * totalQuadOrdinate;
	} else {
	  lk = lk2 + cx * totalQuadOrdinate;
	}
	double piece = lk[qx] * oprob[ox];
	for (int dx=0; dx < maxDims; dx++) piece *= area[dx];
	expected_orig[ox * outRows + cx] += piece;
      }
    }
  }

  int ob_rows;
  int ob_cols;
  getMatrixDims(r_observed_orig, &ob_rows, &ob_cols);
  if (ob_rows != outRows) {
    PrintValue(r_observed_orig);
    pda(expected_orig, outRows, curOutcomes);
    error("Mismatch between observed rows %d and expected rows %d", ob_rows, outRows);
  }
  if (ob_cols != curOutcomes) error("Mismatch between observed cols %d and expected cols %d", ob_cols, curOutcomes);
  int *observed_orig = INTEGER(r_observed_orig);

  for (int rx=0; rx < outRows; rx++) {
    int nk = 0;
    for (int ox=0; ox < curOutcomes; ox++) {
      nk += observed_orig[ox * outRows + rx];
    }
    for (int ox=0; ox < curOutcomes; ox++) {
      expected_orig[ox * outRows + rx] *= nk/denom[rx];
    }
  }

  Free(lk1);
  Free(lk2);

  int df = outRows * (curOutcomes-1);

  const int returnCount = 3;
  SEXP names, ans;
  PROTECT(names = allocVector(STRSXP, returnCount));
  PROTECT(ans = allocVector(VECSXP, returnCount));

  int ansC = -1;
  SET_STRING_ELT(names, ++ansC, mkChar("observed"));
  SET_VECTOR_ELT(ans,   ansC, r_observed_orig);
  SET_STRING_ELT(names, ++ansC, mkChar("expected"));
  SET_VECTOR_ELT(ans,   ansC, r_expected_orig);
  SET_STRING_ELT(names, ++ansC, mkChar("df"));
  SET_VECTOR_ELT(ans,   ansC, ScalarInteger(df));
  if (ansC != returnCount-1) error("Memory corruption");

  namesgets(ans, names);
  UNPROTECT(3);

  return ans;
}

SEXP sumscore_observed(SEXP r_high, SEXP r_data, SEXP r_interest, SEXP r_outcomes)
{
  if (!isInteger(r_data)) error("Data must be of integer type");

  int data_rows;
  int data_cols;
  getMatrixDims(r_data, &data_rows, &data_cols);

  int low = data_cols-1;
  if (low < 1) error("At least 2 columns of data are required");
  int interest = asInteger(r_interest);
  int outcomes = asInteger(r_outcomes);
  int high = asInteger(r_high) - outcomes;

  if (interest < 1 || interest > data_cols)
    error("Interest %d must be between 1 and %d", interest, data_cols);

  int rows = high-low+1;
  SEXP r_ans;
  PROTECT(r_ans = allocMatrix(INTSXP, rows, outcomes));
  int *ans = INTEGER(r_ans);
  memset(ans, 0, sizeof(int) * rows * outcomes);

  int *iresp = INTEGER(r_data) + (interest-1) * data_rows;

  for (int rx=0; rx < data_rows; rx++) {
    int sum=0;
    for (int cx=0; cx < data_cols; cx++) {
      if (cx+1 == interest) continue;
      int *resp = INTEGER(r_data) + cx * data_rows;
      sum += resp[rx];
    }
    ans[(iresp[rx]-1) * rows + sum - low] += 1;
  }

  UNPROTECT(1);
  return r_ans;
}

static double table_concordance(double *mat, int rows, int cols, int ii, int jj)
{
  double sum=0;
  for (int hh=ii+1; hh < rows; ++hh) {
    for (int kk=jj+1; kk < cols; ++kk) {
      sum += mat[kk * rows + hh];
    }
  }
  return sum;
}

static double table_discordance(double *mat, int rows, int cols, int ii, int jj)
{
  double sum=0;
  for (int hh=ii+1; hh < rows; ++hh) {
    for (int kk=0; kk < jj; ++kk) {
      sum += mat[kk * rows + hh];
    }
  }
  return sum;
}

/* See Agresti (1990, p. 22) */
SEXP gamma_cor(SEXP r_mat)
{
  int rows;
  int cols;
  getMatrixDims(r_mat, &rows, &cols);
  SEXP realmat;
  PROTECT(realmat = coerceVector(r_mat, REALSXP));
  double *mat = REAL(realmat);

  double concord = 0;
  for (int ii=0; ii < rows-1; ++ii) {
    for (int jj=0; jj < cols-1; ++jj) {
      concord += mat[jj * rows + ii] * table_concordance(mat, rows, cols, ii, jj);
    }
  }

  double discord = 0;
  for (int ii=0; ii < rows-1; ++ii) {
    for (int jj=1; jj < cols; ++jj) {
      discord += mat[jj * rows + ii] * table_discordance(mat, rows, cols, ii, jj);
    }
  }

  UNPROTECT(1);

  double gamma = (concord - discord) / (concord + discord);
  return ScalarReal(gamma);
}

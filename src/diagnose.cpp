#include "rpf.h"

#include <stdlib.h>

struct ifaGroup {
	bool twotier;
	std::vector<double*> spec;
	int numItems;
	int maxParam;
	double *param;
	int numSpecific;
	int maxAbilities;
	double *mean;
	double *cov;
	std::vector<const char*> factorNames;

	// TODO:
	// scores
	// data

	ifaGroup(bool _twotier) : twotier(_twotier) {}
	void import(SEXP Rlist);
	double *getItemParam(int ix) { return param + maxParam * ix; }
};

void ifaGroup::import(SEXP Rlist)
{
	SEXP argNames;
	Rf_protect(argNames = Rf_getAttrib(Rlist, R_NamesSymbol));

	int mlen = 0;
	int nrow, ncol;
	for (int ax=0; ax < Rf_length(Rlist); ++ax) {
		const char *key = R_CHAR(STRING_ELT(argNames, ax));
		SEXP slotValue = VECTOR_ELT(Rlist, ax);
		if (strEQ(key, "spec")) {
			for (int sx=0; sx < Rf_length(slotValue); ++sx) {
				SEXP model = VECTOR_ELT(slotValue, sx);
				if (!OBJECT(model)) {
					Rf_error("Item models must inherit rpf.base");
				}
				SEXP Rspec;
				Rf_protect(Rspec = R_do_slot(model, Rf_install("spec")));
				spec.push_back(REAL(Rspec));
			}
		} else if (strEQ(key, "param")) {
			param = REAL(slotValue);
			getMatrixDims(slotValue, &maxParam, &numItems);
		} else if (strEQ(key, "mean")) {
			mlen = Rf_length(slotValue);
			mean = REAL(slotValue);
		} else if (strEQ(key, "cov")) {
			getMatrixDims(slotValue, &nrow, &ncol);
			if (nrow != ncol) Rf_error("cov must be a square matrix (not %dx%d)", nrow, ncol);
			cov = REAL(slotValue);

			SEXP dimnames;
			Rf_protect(dimnames = Rf_getAttrib(slotValue, R_DimNamesSymbol));
			if (!Rf_isNull(dimnames) && Rf_length(dimnames) == 2) {
				SEXP names;
				Rf_protect(names = VECTOR_ELT(dimnames, 0));
				int nlen = Rf_length(names);
				factorNames.resize(nlen);
				for (int nx=0; nx < nlen; ++nx) {
					factorNames[nx] = CHAR(STRING_ELT(names, nx));
				}
			}
		} else {
			// ignore
		}
	}
	if (mlen != nrow) Rf_error("Mean length %d does not match cov size %d", mlen, nrow);
	if (numItems != int(spec.size())) {
		Rf_error("param implies %d items but spec is length %d",
			 numItems, (int) spec.size());
	}

	// detect two-tier covariance structure
	std::vector<int> orthogonal;
	if (twotier && mlen >= 3) {
		Eigen::Map<Eigen::MatrixXd> Ecov(cov, mlen, mlen);
		Eigen::Matrix<Eigen::DenseIndex, Eigen::Dynamic, 1> numCov((Ecov.array() != 0.0).matrix().colwise().count());
		std::vector<int> candidate;
		for (int fx=0; fx < numCov.rows(); ++fx) {
			if (numCov(fx) == 1) candidate.push_back(fx);
		}
		if (candidate.size() > 1) {
			std::vector<bool> mask(numItems);
			for (int cx=candidate.size() - 1; cx >= 0; --cx) {
				std::vector<bool> loading(numItems);
				for (int ix=0; ix < numItems; ++ix) {
					loading[ix] = param[ix * maxParam + candidate[cx]] != 0;
				}
				std::vector<bool> overlap(loading.size());
				std::transform(loading.begin(), loading.end(),
					       mask.begin(), overlap.begin(),
					       std::logical_and<bool>());
				if (std::find(overlap.begin(), overlap.end(), true) == overlap.end()) {
					std::transform(loading.begin(), loading.end(),
						       mask.begin(), mask.begin(),
						       std::logical_or<bool>());
					orthogonal.push_back(candidate[cx]);
				}
			}
		}
	}
	std::reverse(orthogonal.begin(), orthogonal.end());
	if (orthogonal.size() && orthogonal[0] != mlen - int(orthogonal.size())) {
		Rf_error("Independent factors must be given after dense factors");
	}

	maxAbilities = mlen;
	numSpecific = orthogonal.size();
}

class ssEAP {
	int lastItem;
public:
	ifaGroup grp;
	ba81NormalQuad quad;
	int maxScore;
	std::vector<int> items;
	std::vector<int> itemOutcomes;
	Eigen::MatrixXd slCur;
	Eigen::MatrixXd slPrev;
	
	ssEAP() : grp(false) {}
	void setup(SEXP grp, double qwidth, int qpts);
	void setLastItem(int which);
	void tpbw1995();
};

#include <iostream> //TODO remove

void ssEAP::setup(SEXP robj, double qwidth, int qpts)
{
	lastItem = 0;

	grp.import(robj);

	Eigen::Map<Eigen::MatrixXd> fullCov(grp.cov, grp.maxAbilities, grp.maxAbilities);
	int dense = grp.maxAbilities - grp.numSpecific;
	Eigen::MatrixXd priCov = fullCov.block(0, 0, dense, dense);
	Eigen::VectorXd sVar = fullCov.diagonal().tail(grp.numSpecific);
	quad.setup(qwidth, qpts, grp.mean, priCov, sVar);

	maxScore = 0;
	itemOutcomes.reserve(grp.numItems);
	for (int cx = 0; cx < grp.numItems; cx++) {
		const double *spec = grp.spec[cx];
		int no = spec[RPF_ISpecOutcomes];
		itemOutcomes.push_back(no);
		maxScore += no - 1;
	}
}

void ssEAP::setLastItem(int which)
{
	lastItem = which;
}

void ssEAP::tpbw1995()
{
	items.reserve(grp.numItems);
	for (int ix=0; ix < grp.numItems; ++ix) if (ix != lastItem) items.push_back(ix);
	items.push_back(lastItem);

	slCur.resize(quad.totalQuadPoints, 1+maxScore);
	int curMax = 0;
	int item0 = items[0];
	{
		const double *spec = grp.spec[item0];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(item0);
		int outcomes = itemOutcomes[item0];
		Eigen::VectorXd oprob(outcomes);
		for (int qx=0; qx < quad.totalQuadPoints; ++qx) {
			double *where = quad.wherePrep.data() + qx * quad.maxDims;
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, quad.maxDims-1)];
			}
			librpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int ox=0; ox < outcomes; ++ox) {
				slCur(qx, ox) = oprob[ox];
			}
		}
		curMax += outcomes - 1;
	}

	slPrev.resize(quad.totalQuadPoints, 1+maxScore);

	for (int curItem=1; curItem < int(items.size()); ++curItem) {
		int ix = items[curItem];
		slCur.swap(slPrev);
		const double *spec = grp.spec[ix];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(ix);
		int outcomes = itemOutcomes[ix];
		Eigen::VectorXd oprob(outcomes);
		slCur.topLeftCorner(slCur.rows(), curMax + outcomes).setZero();
		for (int qx=0; qx < quad.totalQuadPoints; ++qx) {
			double *where = quad.wherePrep.data() + qx * quad.maxDims;
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, quad.maxDims-1)];
			}
			librpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int cx=0; cx <= curMax; cx++) {
				for (int ox=0; ox < outcomes; ox++) {
					slCur(qx, cx + ox) += slPrev(qx, cx) * oprob[ox];
				}
			}
		}
		curMax += outcomes - 1;
	}
}

SEXP ot2000_wrapper(SEXP robj, SEXP Ritem, SEXP Rwidth, SEXP Rpts, SEXP Ralter)
{
	omxManageProtectInsanity mpi;

	double qwidth = Rf_asReal(Rwidth);
	int qpts = Rf_asInteger(Rpts);
	int interest = Rf_asInteger(Ritem) - 1;

	ssEAP myeap;
	myeap.setup(robj, qwidth, qpts);
	myeap.setLastItem(interest);
	myeap.tpbw1995();

	int outcomes = myeap.itemOutcomes[interest];

	ifaGroup &grp = myeap.grp;
	ba81NormalQuad &quad = myeap.quad;
	Eigen::MatrixXd iProb(quad.totalQuadPoints, outcomes);

	{
		const double *spec = grp.spec[interest];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(interest);
		Eigen::VectorXd oprob(outcomes);
		for (int qx=0; qx < quad.totalQuadPoints; ++qx) {
			double *where = quad.wherePrep.data() + qx * quad.maxDims;
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, quad.maxDims-1)];
			}
			librpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int ox=0; ox < outcomes; ox++) {
				iProb(qx, ox) = oprob[ox];
			}
		}
	}

	if (Rf_asLogical(Ralter)) {
		Eigen::MatrixXd &slCur = myeap.slCur;
		Eigen::MatrixXd &slPrev = myeap.slPrev;

		Eigen::Map<Eigen::VectorXd> areaCol(quad.priQarea.data(), quad.priQarea.size());
		slCur.array().colwise() *= areaCol.array();
		slPrev.array().colwise() *= areaCol.array();

		Eigen::VectorXd ssProb(1+myeap.maxScore);
		ssProb = slCur.array().colwise().sum();

		SEXP Rexpected;
		Rf_protect(Rexpected = Rf_allocMatrix(REALSXP, 1+myeap.maxScore, outcomes));
		double *out = REAL(Rexpected);
		memset(out, 0, sizeof(double) * (1+myeap.maxScore) * outcomes);

		for (int startScore=0; startScore+outcomes <= 1+myeap.maxScore; ++startScore) {
			Eigen::MatrixXd numerDen(quad.totalQuadPoints, outcomes);
			numerDen = iProb.array().colwise() * slPrev.col(startScore).array();
			Eigen::VectorXd numer(outcomes);
			numer = numerDen.colwise().sum();
			for (int ox=0; ox < outcomes; ++ox) {
				out[ox * (1+myeap.maxScore) + (startScore+ox)] = numer(ox) / ssProb(startScore + ox);
			}
		}

		return Rexpected;
	} else {
		Eigen::MatrixXd &slPrev = myeap.slPrev;
	
		Eigen::Map<Eigen::VectorXd> areaCol(quad.priQarea.data(), quad.priQarea.size());
		slPrev.array().colwise() *= areaCol.array();

		Eigen::VectorXd ssProb(1+myeap.maxScore);
		ssProb = slPrev.array().colwise().sum();

		int prevMaxScore = myeap.maxScore - (outcomes-1);
		SEXP Rexpected;
		Rf_protect(Rexpected = Rf_allocMatrix(REALSXP, prevMaxScore + 1, outcomes));
		double *out = REAL(Rexpected);

		for (int startScore=0; startScore <= prevMaxScore; ++startScore) {
			Eigen::MatrixXd numerDen(quad.totalQuadPoints, outcomes);
			numerDen = iProb.array().colwise() * slPrev.col(startScore).array();
			Eigen::VectorXd numer(outcomes);
			numer = numerDen.colwise().sum();
			for (int ox=0; ox < outcomes; ++ox) {
				out[ox * (prevMaxScore+1) + startScore] = numer(ox) / ssProb(startScore);
			}
		}

		return Rexpected;
	}
}

SEXP sumscoreEAP(SEXP robj, SEXP Rwidth, SEXP Rpts)
{
	omxManageProtectInsanity mpi;

	double qwidth = Rf_asReal(Rwidth);
	int qpts = Rf_asInteger(Rpts);

	ssEAP myeap;
	myeap.setup(robj, qwidth, qpts);
	myeap.tpbw1995();

	ba81NormalQuad &quad = myeap.quad;
	ifaGroup &grp = myeap.grp;
	Eigen::MatrixXd &slCur = myeap.slCur;

	if (quad.numSpecific == 0) {
		Eigen::Map<Eigen::VectorXd> areaCol(quad.priQarea.data(), quad.priQarea.size());
		slCur.array().colwise() *= areaCol.array();
	} else {
		Eigen::VectorXd den(quad.totalQuadPoints);
		for (int qloc=0, qx=0; qx < quad.totalPrimaryPoints; ++qx) {
			for (int sx=0; sx < quad.quadGridSize; ++sx) {
				den(qloc) = quad.priQarea[qx] * quad.speQarea[sx * quad.numSpecific];
				++qloc;
			}
		}
		slCur.array().colwise() *= den.array();
	}

	int curMax = myeap.maxScore;
	int outRows = 1 + curMax;
	int outCols = 1 + 2 * quad.maxAbilities + triangleLoc1(quad.maxAbilities);

	SEXP dimnames;
	Rf_protect(dimnames = Rf_allocVector(VECSXP, 2));
	SEXP names;
	Rf_protect(names = Rf_allocVector(STRSXP, outRows));
	for (int rx=0; rx <= curMax; ++rx) {
		const int SMALLBUF = 10;
		char buf[SMALLBUF];
		snprintf(buf, SMALLBUF, "%d", rx);
                SET_STRING_ELT(names, rx, Rf_mkChar(buf));
	}
	SET_VECTOR_ELT(dimnames, 0, names);

	Rf_protect(names = Rf_allocVector(STRSXP, outCols));
	SET_STRING_ELT(names, 0, Rf_mkChar("p"));
	for (int ax=0; ax < quad.maxAbilities; ++ax) {
		const int SMALLBUF = 10;
		char buf[SMALLBUF];
		if (grp.factorNames.size()) {
			SET_STRING_ELT(names, 1+ax, Rf_mkChar(grp.factorNames[ax]));
		} else {
			snprintf(buf, SMALLBUF, "f%d", 1+ax);
			SET_STRING_ELT(names, 1+ax, Rf_mkChar(buf));
		}
		snprintf(buf, SMALLBUF, "se%d", 1+ax);
		SET_STRING_ELT(names, 1+quad.maxAbilities+ax, Rf_mkChar(buf));
	}
	for (int cx=0; cx < triangleLoc1(quad.maxAbilities); ++cx) {
		const int SMALLBUF = 10;
		char buf[SMALLBUF];
		snprintf(buf, SMALLBUF, "cov%d", 1+cx);
		SET_STRING_ELT(names, 1+2*quad.maxAbilities+cx, Rf_mkChar(buf));
	}
	SET_VECTOR_ELT(dimnames, 1, names);

	Eigen::VectorXd ssProb(outRows);
	ssProb = slCur.array().colwise().sum();

	SEXP Rout;
	Rf_protect(Rout = Rf_allocMatrix(REALSXP, outRows, outCols));
	Rf_setAttrib(Rout, R_DimNamesSymbol, dimnames);
	double *out = REAL(Rout);
	memcpy(out, ssProb.data(), sizeof(double) * outRows);
	for (int cx=0; cx <= curMax; cx++) {
		std::vector<double> pad(quad.maxAbilities + triangleLoc1(quad.maxAbilities));
		if (quad.numSpecific == 0) {
			quad.EAP(&slCur.coeffRef(0, cx), 1/ssProb[cx], pad.data());
		} else {
			Eigen::MatrixXd weight(quad.numSpecific, quad.totalQuadPoints);
			weight.array().rowwise() = slCur.array().col(cx).transpose();
			quad.EAP(weight.data(), 1/ssProb[cx], pad.data());
		}
		for (int sx=0; sx < quad.maxAbilities; ++sx) {
			out[(1+sx) * outRows + cx] = pad[sx];
			out[(1+quad.maxAbilities+sx) * outRows + cx] = sqrt(pad[quad.maxAbilities + triangleLoc0(sx)]);
		}
		for (int sx=0; sx < triangleLoc1(quad.maxAbilities); ++sx) {
			out[(1+2*quad.maxAbilities+sx) * outRows + cx] = pad[quad.maxAbilities + sx];
		}
	}

	return Rout;
}

SEXP pairwiseExpected(SEXP robj, SEXP Rwidth, SEXP Rpts, SEXP Ritems)
{
	omxManageProtectInsanity mpi;

	double qwidth = Rf_asReal(Rwidth);
	int qpts = Rf_asInteger(Rpts);
	if (Rf_length(Ritems) != 2) Rf_error("A pair of items must be specified");

	ifaGroup grp(true);
	grp.import(robj);
	
	// factor out? TODO
	Eigen::Map<Eigen::MatrixXd> fullCov(grp.cov, grp.maxAbilities, grp.maxAbilities);
	int dense = grp.maxAbilities - grp.numSpecific;
	Eigen::MatrixXd priCov = fullCov.block(0, 0, dense, dense);
	Eigen::VectorXd sVar = fullCov.diagonal().tail(grp.numSpecific);

	ba81NormalQuad quad;
	quad.setup(qwidth, qpts, grp.mean, priCov, sVar);

	int i1 = INTEGER(Ritems)[0];
	int i2 = INTEGER(Ritems)[1];
	if (i1 < 0 || i1 >= (int) grp.spec.size()) Rf_error("Item %d out of range", i1);
	if (i2 < 0 || i2 >= (int) grp.spec.size()) Rf_error("Item %d out of range", i2);
	if (i1 == i2) Rf_warning("Request to create bivariate distribution of %d with itself", i1);

	double *i1par = &grp.param[i1 * grp.maxParam];
	double *i2par = &grp.param[i2 * grp.maxParam];

	int specific = -1;
	if (grp.numSpecific) {
		for (int ax=quad.maxDims; ax < quad.maxAbilities; ax++) {
			if (i1par[ax] != 0 && i2par[ax] != 0) {
				specific = ax - quad.maxDims;
				break;
			}
		}
	}

	double *spec1 = grp.spec[i1];
	int id1 = spec1[RPF_ISpecID];
	int outcomes1 = spec1[RPF_ISpecOutcomes];

	double *spec2 = grp.spec[i2];
	int id2 = spec1[RPF_ISpecID];
	int outcomes2 = spec2[RPF_ISpecOutcomes];

	SEXP Rexpected;
	Rf_protect(Rexpected = Rf_allocMatrix(REALSXP, outcomes1, outcomes2));
	double *outMem = REAL(Rexpected);
	Eigen::Map<Eigen::MatrixXd> out(outMem, outcomes1, outcomes2);
	out.setZero();

	Eigen::VectorXd o1(outcomes1);
	Eigen::VectorXd o2(outcomes2);

	if (specific == -1) {
		for (int qx=0; qx < quad.totalQuadPoints; ++qx) {
			double *where = quad.wherePrep.data() + qx * quad.maxDims;
			(*librpf_model[id1].prob)(spec1, i1par, where, o1.data());
			(*librpf_model[id2].prob)(spec2, i2par, where, o2.data());
			out += (o1 * o2.transpose()) * quad.priQarea[qx];
		}
	} else {
		Eigen::VectorXd ptheta(quad.maxAbilities);
		for (int qloc=0, qx=0; qx < quad.totalPrimaryPoints; ++qx) {
			for (int sx=0; sx < quad.quadGridSize; ++sx) {
				double *where = quad.wherePrep.data() + qloc * quad.maxDims;
				for (int dx=0; dx < quad.maxAbilities; dx++) {
					ptheta[dx] = where[std::min(dx, quad.maxDims-1)];
				}
				(*librpf_model[id1].prob)(spec1, i1par, ptheta.data(), o1.data());
				(*librpf_model[id2].prob)(spec2, i2par, ptheta.data(), o2.data());
				double area = quad.priQarea[qx] * quad.speQarea[sx * quad.numSpecific + specific];
				out += (o1 * o2.transpose()) * area;
				++qloc;
			}
		}
	}

	return Rexpected;
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

struct ManhattenCollapse {
	Eigen::Map<Eigen::ArrayXXi> obs;
	Eigen::Map<Eigen::ArrayXXd> expected;

	Eigen::DenseIndex smr, smc;
	double bestFit;
	Eigen::DenseIndex bestR, bestC;
	
	ManhattenCollapse(int rows, int cols, int *oMem, double *eMem)
		: obs(oMem, rows, cols), expected(eMem, rows, cols) {}
	void probe(Eigen::DenseIndex br, Eigen::DenseIndex bc);
	int run();
};

void ManhattenCollapse::probe(Eigen::DenseIndex br, Eigen::DenseIndex bc)
{
	if (br < 0 || bc < 0 || br >= expected.rows() || bc >= expected.cols()) return;
	if (expected(br, bc) < bestFit) {
		bestFit = expected(br, bc);
		bestR = br;
		bestC = bc;
	}
}

int ManhattenCollapse::run()
{
	const double worstFit = 1e100;
	const int maxDist = obs.rows() + obs.cols();
	int collapsed = 0;

	while (expected.minCoeff(&smr, &smc) < KANG_CHEN_MIN_EXPECTED) {
		bool done = false;
		for (int dist=1; dist < maxDist && !done; ++dist) {
			for (int updown=0; updown <= dist && !done; ++updown) {
				int leftright = dist - updown;
				bestFit = worstFit;
				probe(smr + updown, smc + leftright);
				probe(smr + updown, smc - leftright);
				probe(smr - updown, smc + leftright);
				probe(smr - updown, smc - leftright);
				if (bestFit < worstFit) {
					expected(bestR, bestC) += expected(smr, smc);
					obs(bestR, bestC) += obs(smr, smc);
					expected(smr, smc) = NA_REAL;
					obs(smr, smc) = NA_REAL;
					++collapsed;
					done = true;
				}
			}
		}
		if (!done) {
			pda(expected.data(), expected.rows(), expected.cols());
			Rf_error("Collapse algorithm failed");
		}
	}

	return collapsed;
}

SEXP kang_chen_2007_wrapper(SEXP r_observed_orig, SEXP r_expected_orig)
{
	omxManageProtectInsanity mpi;

  int rows, cols;
  getMatrixDims(r_expected_orig, &rows, &cols);

  {
    int orows, ocols;
    getMatrixDims(r_observed_orig, &orows, &ocols);
    if (rows != orows || cols != ocols)
	    Rf_error("Observed %dx%d and expected %dx%d matrices must have same dimensions",
		     orows, ocols, rows, cols);
  }

  SEXP r_observed, r_expected;
  Rf_protect(r_observed = Rf_allocMatrix(INTSXP, rows, cols));
  Rf_protect(r_expected = Rf_allocMatrix(REALSXP, rows, cols));

  int *observed = INTEGER(r_observed);
  double *expected = REAL(r_expected);
  memcpy(observed, INTEGER(r_observed_orig), sizeof(int) * rows * cols);
  memcpy(expected, REAL(r_expected_orig), sizeof(double) * rows * cols);

  ManhattenCollapse mcollapse(rows, cols, observed, expected);
  int collapsed = mcollapse.run();

  MxRList out;
  out.add("O", r_observed);
  out.add("E", r_expected);
  out.add("collapsed", Rf_ScalarInteger(collapsed));
  return out.asR();
}

SEXP sumscore_observed(SEXP r_high, SEXP r_data, SEXP r_interest, SEXP r_outcomes, SEXP Ralter)
{
	bool alter = Rf_asLogical(Ralter);
  int data_rows = Rf_length(VECTOR_ELT(r_data, 0));
  int data_cols = Rf_length(r_data);

  int interest = Rf_asInteger(r_interest);
  int outcomes = Rf_asInteger(r_outcomes);
  int high;
  if (!alter) {
	  high = Rf_asInteger(r_high) - (outcomes - 1);
  } else {
	  high = Rf_asInteger(r_high);
  }

  if (interest < 1 || interest > data_cols)
    Rf_error("Interest %d must be between 1 and %d", interest, data_cols);

  SEXP r_ans;
  Rf_protect(r_ans = Rf_allocMatrix(INTSXP, high, outcomes));
  int *ans = INTEGER(r_ans);
  memset(ans, 0, sizeof(int) * high * outcomes);

  int *iresp = INTEGER(VECTOR_ELT(r_data, interest-1));

  for (int rx=0; rx < data_rows; rx++) {
	  bool missing = false;
    int sum=0;
    for (int cx=0; cx < data_cols; cx++) {
      if (!alter && cx+1 == interest) continue;
      int *resp = INTEGER(VECTOR_ELT(r_data, cx));
      if (resp[rx] == NA_INTEGER) {
	      missing = true;
	      break;
      }
      sum += resp[rx] - 1;
    }
    if (missing) continue;
    int pick = iresp[rx];
    if (pick == NA_INTEGER) continue;
    ans[(pick-1) * high + sum] += 1;
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
  Rf_protect(realmat = Rf_coerceVector(r_mat, REALSXP));
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
  return Rf_ScalarReal(gamma);
}

template <typename T1, typename T2, typename T3>
static inline double crosstabMS(Eigen::ArrayBase<T1> &observed,
				Eigen::ArrayBase<T2> &expected,
				Eigen::ArrayBase<T3> &rowSum)
{
	Eigen::ArrayXXd diff(observed.rows(), observed.cols());
	observed.colwise() /= rowSum;
	diff = observed - expected;
	if (observed.rows() == 1) {
		return ((diff * diff).rowwise().sum() * rowSum).sum();
	} else {
		// not sure if this is correct TODO
		diff.colwise() *= rowSum;
		return (diff * diff).sum();
	}
}

SEXP crosstabTest(SEXP Robserved, SEXP Rexpected, SEXP Rtrials)
{
	int rows, cols;
	getMatrixDims(Robserved, &rows, &cols);
	{
		int erows, ecols;
		getMatrixDims(Rexpected, &erows, &ecols);
		if (rows != erows || cols != ecols) {
			Rf_error("observed and expected matrices must be the same dimension");
		}
	}

	Eigen::ArrayXXd observed(rows, cols);
	if (Rf_isInteger(Robserved)) {
		Eigen::Map<Eigen::ArrayXXi > tmp(INTEGER(Robserved), rows, cols);
		observed = tmp.cast<double>();
	} else if (Rf_isReal(Robserved)) {
		memcpy(observed.data(), REAL(Robserved), sizeof(double) * rows * cols);
	} else {
		Rf_error("observed is an unknown type");
	}

	Eigen::ArrayXXd expected(rows, cols);
	memcpy(expected.data(), REAL(Rexpected), sizeof(double) * rows * cols);

	Eigen::ArrayXd rowSum(rows);
	rowSum = observed.rowwise().sum();
	if (((expected.rowwise().sum() - rowSum).abs() > 1e-6).any()) {
		Rf_error("observed and expected row sums must match");
	}

	expected.colwise() /= rowSum;

	Eigen::VectorXi simSize(rows);
	simSize = rowSum.cast<int>();

	if (rows == 1) {
		// Perkins, Tygert, Ward (2011, p. 12)
		for (int sx=0; sx < rows; ++sx) simSize(sx) = std::min(simSize(sx), 185);
	}
	Eigen::ArrayXd simSizeD(rows);
	simSizeD = simSize.cast<double>();

	double refMS = crosstabMS(observed, expected, rowSum);

	Eigen::ArrayXXi Eprob(rows, cols);
	Eprob = (expected * RAND_MAX).cast<int>();

	int trials = Rf_asInteger(Rtrials); // SE = 1/(.05 * .95 * sqrt(trials))
	Eigen::ArrayXd mcMS(trials);
	for (int tx=0; tx < trials; ++tx) {
		Eigen::ArrayXXd draw(rows, cols);
		draw.setZero();
		for (int rx=0; rx < rows; ++rx) {
			for (int sx=0; sx < simSize(rx); ++sx) {
				int r1 = rand();
				for (int cx=0; cx < cols; ++cx) {
					int threshold = Eprob(rx, cx);
					if (r1 <= threshold) {
						draw(rx, cx) += 1;
						break;
					} else {
						r1 -= threshold;
					}
				}
			}
		}
		mcMS(tx) = crosstabMS(draw, expected, simSizeD);
	}
	Eigen::DenseIndex cnt = (mcMS >= refMS).count();
	return Rf_ScalarReal(double(cnt) / trials);
}

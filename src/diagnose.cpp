#include "rpf.h"

class ssEAP {
	int lastItem;
	void tpbw1995Prep();

	template <typename T1, typename T2>
	void aggregateSpecific(Eigen::ArrayBase<T1> &inMat, Eigen::ArrayBase<T2> &Eis);

	template <typename T1, typename T2, typename T3>
	void tt2prod(Eigen::ArrayBase<T1> &slCur, Eigen::ArrayBase<T2> &buffer, Eigen::ArrayBase<T3> &perSpecific);

	template <typename T1, typename T2, typename T3>
	void prod2ss(Eigen::ArrayBase<T1> &buffer, Eigen::ArrayBase<T2> &ssCur, Eigen::ArrayBase<T3> &perSpecific);

public:
	ifaGroup grp;
	int *mask;
	int maxScore;
	std::vector<int> items;

	Eigen::ArrayXXd ttCur;
	Eigen::ArrayXi ttCurMax;
	Eigen::ArrayXXd slCur;
	Eigen::ArrayXd ssProbCur;

	Eigen::ArrayXXd ttPrev;
	Eigen::ArrayXi ttPrevCurMax;
	Eigen::ArrayXXd slPrev;
	Eigen::ArrayXd ssProbPrev;
	
	ssEAP(bool twotier) : grp(1, twotier) {}
	void setup(SEXP grp, double qwidth, int qpts, int *_mask);
	void setLastItem(int which);
	void tpbw1995Vanilla();
	void tpbw1995TwoTier();
	void tpbw1995();

	// tt2ss == two-tier to sum-score
	template <typename T1, typename T2, typename T3>
	void tt2ss(Eigen::ArrayBase<T1> &curMax, Eigen::ArrayBase<T2> &slCur, Eigen::ArrayBase<T3> &ssProbCur);
};

void ssEAP::setup(SEXP robj, double qwidth, int qpts, int *_mask)
{
	lastItem = -1;
	mask = _mask;

	grp.setGridFineness(qwidth, qpts);
	grp.import(robj, true);
	grp.setupQuadrature();
}

void ssEAP::setLastItem(int which)
{
	lastItem = which;
}

void ssEAP::tpbw1995Prep()
{
	maxScore = 0;
	for (int cx = 0; cx < grp.numItems(); cx++) {
		const double *spec = grp.spec[cx];
		int no = spec[RPF_ISpecOutcomes];
		if ((lastItem != -1 && cx == lastItem) || mask[cx]) {
			maxScore += no - 1;
			if (cx != lastItem) items.push_back(cx);
		}
	}

	if (lastItem >= 0) items.push_back(lastItem);
}

void ssEAP::tpbw1995()
{
	tpbw1995Prep();

	if (grp.numSpecific == 0) {
		tpbw1995Vanilla();
	} else {
		tpbw1995TwoTier();
	}
}

void ssEAP::tpbw1995Vanilla()
{
	ba81NormalQuad &quad = grp.quad;
	slCur.resize(quad.totalQuadPoints, 1+maxScore);
	Eigen::Map<Eigen::ArrayXd> areaCol(quad.priQarea.data(), quad.priQarea.size());
	slCur.colwise() = areaCol;

	int curMax = 0;
	int item0 = items[0];
	{
		const double *spec = grp.spec[item0];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(item0);
		int outcomes = grp.itemOutcomes[item0];
		Eigen::VectorXd oprob(outcomes);
		for (int qx=0; qx < quad.totalQuadPoints; ++qx) {
			double *where = quad.wherePrep.data() + qx * quad.maxDims;
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, quad.maxDims-1)];
			}
			Glibrpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int ox=0; ox < outcomes; ++ox) {
				slCur(qx, ox) *= oprob[ox];
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
		int outcomes = grp.itemOutcomes[ix];
		Eigen::VectorXd oprob(outcomes);
		slCur.topLeftCorner(slCur.rows(), curMax + outcomes).setZero();
		for (int qx=0; qx < quad.totalQuadPoints; ++qx) {
			double *where = quad.wherePrep.data() + qx * quad.maxDims;
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, quad.maxDims-1)];
			}
			Glibrpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int cx=0; cx <= curMax; cx++) {
				for (int ox=0; ox < outcomes; ox++) {
					slCur(qx, cx + ox) += slPrev(qx, cx) * oprob[ox];
				}
			}
		}
		curMax += outcomes - 1;
		R_CheckUserInterrupt();
	}

	ssProbCur.resize(1+maxScore);
	ssProbCur = slCur.colwise().sum();
	ssProbPrev.resize(1+maxScore);     // make smaller TODO
	ssProbPrev = slPrev.colwise().sum();
}

template <typename T1, typename T2>
void ssEAP::aggregateSpecific(Eigen::ArrayBase<T1> &inMat, Eigen::ArrayBase<T2> &Eis)
{
	ba81NormalQuad &quad = grp.quad;

	Eis.setZero();
	for (int qx=0, qloc = 0, eisloc = 0; qx < quad.totalPrimaryPoints; qx++) {
		for (int sx=0; sx < quad.quadGridSize; sx++) {
			for (int sgroup=0; sgroup < quad.numSpecific; ++sgroup) {
				Eis.row(eisloc + sgroup) += inMat.row(qloc);
				++qloc;
			}
		}
		eisloc += quad.numSpecific;
	}
}

template <typename T1, typename T2, typename T3>
void ssEAP::tt2prod(Eigen::ArrayBase<T1> &tt, Eigen::ArrayBase<T2> &buffer, Eigen::ArrayBase<T3> &perSpecific)
{
	ba81NormalQuad &quad = grp.quad;

	int combinations = perSpecific.prod();
	int destRows = tt.rows() / perSpecific.rows();
	buffer.setOnes();

	for (int qx=0; qx < destRows; qx++) {
		for (int cx=0; cx < combinations; ++cx) {
			int chip = cx;
			for (int sgroup=0; sgroup < quad.numSpecific; ++sgroup) {
				int col = chip % perSpecific[sgroup];
				chip /= perSpecific[sgroup];
				buffer(qx, cx) *= tt(qx * perSpecific.rows() + sgroup, col);
			}
		}
	}
}

template <typename T1, typename T2, typename T3>
void ssEAP::prod2ss(Eigen::ArrayBase<T1> &buffer, Eigen::ArrayBase<T2> &ssMat, Eigen::ArrayBase<T3> &perSpecific)
{
	ba81NormalQuad &quad = grp.quad;

	int combinations = perSpecific.prod();

	ssMat.setZero();
	for (int cx=0; cx < combinations; ++cx) {
		int chip = cx;
		int ss = 0;
		for (int sgroup=0; sgroup < quad.numSpecific; ++sgroup) {
			ss += chip % perSpecific[sgroup];
			chip /= perSpecific[sgroup];
		}
		ssMat.col(ss) += buffer.col(cx);
	}
}

template <typename T1, typename T2, typename T3>
void ssEAP::tt2ss(Eigen::ArrayBase<T1> &curMax1, Eigen::ArrayBase<T2> &curTbl,
		  Eigen::ArrayBase<T3> &outTbl)
{
	int numScores = 1+maxScore;
	ba81NormalQuad &quad = grp.quad;

	Eigen::ArrayXXd Eis(quad.totalPrimaryPoints * quad.numSpecific, numScores);
	aggregateSpecific(curTbl, Eis);

	Eigen::ArrayXi perSpecific = curMax1 + 1;
	int combinations = perSpecific.prod();

	Eigen::ArrayXXd prodEis(quad.totalPrimaryPoints, combinations);
	tt2prod(Eis, prodEis, perSpecific);

	outTbl.derived().resize(quad.totalPrimaryPoints, numScores);
	prod2ss(prodEis, outTbl, perSpecific);

	Eigen::Map<Eigen::ArrayXd> areaCol(quad.priQarea.data(), quad.priQarea.size());
	outTbl.colwise() *= areaCol;
}

void ssEAP::tpbw1995TwoTier()
{
	int numScores = 1+maxScore;

	ba81NormalQuad &quad = grp.quad;

	Eigen::ArrayXd speAreaCol(quad.totalQuadPoints * quad.numSpecific);
	for (int qx=0, qloc = 0; qx < quad.totalPrimaryPoints; qx++) {
		for (int sx=0; sx < quad.quadGridSize * quad.numSpecific; sx++) {
			speAreaCol[qloc] = quad.speQarea[sx];
			++qloc;
		}
	}

	ttCur.resize(quad.totalQuadPoints * quad.numSpecific, numScores);
	ttPrev.resize(quad.totalQuadPoints * quad.numSpecific, numScores);
	ttCur.colwise() = speAreaCol;
	ttPrev.colwise() = speAreaCol;

	ttCurMax.resize(quad.numSpecific);
	ttCurMax.setZero();
	ttPrevCurMax = ttCurMax;

	int item0 = items[0];
	{
		const double *spec = grp.spec[item0];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(item0);
		int outcomes = grp.itemOutcomes[item0];
		int Sgroup = grp.Sgroup[item0];
		Eigen::VectorXd oprob(outcomes);
		for (int qx=0; qx < quad.totalQuadPoints; ++qx) {
			double *where = quad.wherePrep.data() + qx * quad.maxDims;
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, quad.maxDims-1)];
			}
			Glibrpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int ox=0; ox < outcomes; ++ox) {
				ttCur(qx * quad.numSpecific + Sgroup, ox) *= oprob[ox];
			}
		}
		ttCurMax(Sgroup) += outcomes - 1;
	}

	for (int curItem=1; curItem < int(items.size()); ++curItem) {
		int ix = items[curItem];
		ttPrev = ttCur; // can't swap because we only update the item's Sgroup
		ttPrevCurMax = ttCurMax;
		const double *spec = grp.spec[ix];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(ix);
		int outcomes = grp.itemOutcomes[ix];
		int Sgroup = grp.Sgroup[ix];
		Eigen::VectorXd oprob(outcomes);
		for (int qx=0; qx < quad.totalQuadPoints; ++qx) {
			int row = qx * quad.numSpecific + Sgroup;
			ttCur.row(row).setZero();
			double *where = quad.wherePrep.data() + qx * quad.maxDims;
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, quad.maxDims-1)];
			}
			Glibrpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int cx=0; cx <= ttCurMax(Sgroup); cx++) {
				for (int ox=0; ox < outcomes; ox++) {
					ttCur(row, cx + ox) += ttPrev(row, cx) * oprob[ox];
				}
			}
		}
		ttCurMax(Sgroup) += outcomes - 1;
		R_CheckUserInterrupt();
	}

	tt2ss(ttCurMax, ttCur, slCur);
	ssProbCur.resize(numScores);
	ssProbCur = slCur.colwise().sum();

	tt2ss(ttPrevCurMax, ttPrev, slPrev);
	ssProbPrev.resize(numScores);
	ssProbPrev = slPrev.colwise().sum();
}

template <typename T1, typename T2>
void otMix(ssEAP &myeap, int Sgroup, int ox, Eigen::ArrayBase<T1> &iProb, Eigen::ArrayBase<T2> &numer)
{
	ba81NormalQuad &quad = myeap.grp.quad;

	if (quad.numSpecific == 0) {
		numer = (myeap.slPrev.colwise() * iProb.col(ox)).colwise().sum();
	} else {
		Eigen::ArrayXXd ttPrev = myeap.ttPrev;
		for (int qx = 0; qx < quad.totalQuadPoints; ++qx) {
			ttPrev.row(qx*quad.numSpecific + Sgroup) *= iProb(qx, ox);
		}
		Eigen::ArrayXXd ssPrev;
		myeap.tt2ss(myeap.ttPrevCurMax, ttPrev, ssPrev);
		numer = ssPrev.colwise().sum();
	}
}

SEXP ot2000_wrapper(SEXP robj, SEXP Ritem, SEXP Rwidth, SEXP Rpts, SEXP Ralter,
		    SEXP Rmask, SEXP Rtwotier)
{
	ProtectAutoBalanceDoodad mpi;

	double qwidth = Rf_asReal(Rwidth);
	int qpts = Rf_asInteger(Rpts);
	int interest = Rf_asInteger(Ritem) - 1;

	ssEAP myeap(Rf_asLogical(Rtwotier));
	myeap.setup(robj, qwidth, qpts, LOGICAL(Rmask));
	myeap.setLastItem(interest);
	myeap.tpbw1995();

	int outcomes = myeap.grp.itemOutcomes[interest];

	ifaGroup &grp = myeap.grp;
	ba81NormalQuad &quad = myeap.grp.quad;
	Eigen::ArrayXXd iProb(quad.totalQuadPoints, outcomes);

	{
		const double *spec = grp.spec[interest];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(interest);
		Eigen::ArrayXd oprob(outcomes);
		for (int qx=0; qx < quad.totalQuadPoints; ++qx) {
			double *where = quad.wherePrep.data() + qx * quad.maxDims;
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, quad.maxDims-1)];
			}
			Glibrpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int ox=0; ox < outcomes; ox++) {
				iProb(qx, ox) = oprob[ox];
			}
		}
	}

	int Sgroup = 0;
	if (quad.numSpecific) Sgroup = grp.Sgroup[interest];

	if (Rf_asLogical(Ralter)) {
		// as documented in various publications
		Eigen::ArrayXd &ssProb = myeap.ssProbCur;

		SEXP Rexpected;
		Rf_protect(Rexpected = Rf_allocMatrix(REALSXP, 1+myeap.maxScore, outcomes));
		Eigen::Map< Eigen::ArrayXXd > out(REAL(Rexpected), 1+myeap.maxScore, outcomes);
		out.setZero();

		for (int ox=0; ox < outcomes; ++ox) {
			Eigen::ArrayXd numer;
			otMix(myeap, Sgroup, ox, iProb, numer);
			for (int startScore=0; startScore+outcomes <= 1+myeap.maxScore; ++startScore) {
				out(startScore+ox, ox) = numer(startScore) / ssProb(startScore + ox);
			}
		}

		return Rexpected;
	} else {
		// slightly more powerful for small number of items (IRTPRO/flexMIRT)
		Eigen::ArrayXd &ssProb = myeap.ssProbPrev;

		int prevMaxScore = myeap.maxScore - (outcomes-1);
		SEXP Rexpected;
		Rf_protect(Rexpected = Rf_allocMatrix(REALSXP, prevMaxScore + 1, outcomes));
		Eigen::Map< Eigen::ArrayXXd > out(REAL(Rexpected), prevMaxScore + 1, outcomes);

		for (int ox=0; ox < outcomes; ++ox) {
			Eigen::ArrayXd numer;
			otMix(myeap, Sgroup, ox, iProb, numer);
			for (int startScore=0; startScore <= prevMaxScore; ++startScore) {
				out(startScore, ox) = numer(startScore) / ssProb(startScore);
			}
		}

		return Rexpected;
	}
}

SEXP sumscoreEAP(SEXP robj, SEXP Rwidth, SEXP Rpts, SEXP Rmask, SEXP twotier, SEXP debug)
{
	ProtectAutoBalanceDoodad mpi;

	double qwidth = Rf_asReal(Rwidth);
	int qpts = Rf_asInteger(Rpts);

	ssEAP myeap(Rf_asLogical(twotier));
	myeap.setup(robj, qwidth, qpts, LOGICAL(Rmask));
	myeap.tpbw1995();

	ba81NormalQuad &quad = myeap.grp.quad;
	ifaGroup &grp = myeap.grp;
	Eigen::ArrayXXd &slCur = myeap.slCur;

	if (Rf_asLogical(debug)) {
		SEXP Rout;
		Rf_protect(Rout = Rf_allocMatrix(REALSXP, slCur.rows(), slCur.cols()));
		memcpy(REAL(Rout), slCur.data(), sizeof(double) * slCur.rows() * slCur.cols());
		return Rout;
	}

	int curMax = myeap.maxScore;
	int outRows = 1 + curMax;
	int outCols = 1 + 2 * quad.primaryDims + triangleLoc1(quad.primaryDims);

	SEXP dimnames;
	Rf_protect(dimnames = Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dimnames, 0, R_NilValue);

	SEXP names;
	Rf_protect(names = Rf_allocVector(STRSXP, outCols));
	SET_STRING_ELT(names, 0, Rf_mkChar("p"));
	for (int ax=0; ax < quad.primaryDims; ++ax) {
		const int SMALLBUF = 20;
		char buf[SMALLBUF];
		SET_STRING_ELT(names, 1+ax, Rf_mkChar(grp.factorNames[ax].c_str()));
		snprintf(buf, SMALLBUF, "se%d", 1+ax);
		SET_STRING_ELT(names, 1+quad.primaryDims+ax, Rf_mkChar(buf));
	}
	for (int cx=0; cx < triangleLoc1(quad.primaryDims); ++cx) {
		const int SMALLBUF = 20;
		char buf[SMALLBUF];
		snprintf(buf, SMALLBUF, "cov%d", 1+cx);
		SET_STRING_ELT(names, 1+2*quad.primaryDims+cx, Rf_mkChar(buf));
	}
	SET_VECTOR_ELT(dimnames, 1, names);

	Eigen::ArrayXd &ssProb = myeap.ssProbCur;

	ba81NormalQuad pquad;  //primary only
	Eigen::Map<Eigen::MatrixXd> fullCov(grp.cov, grp.maxAbilities, grp.maxAbilities);
	Eigen::MatrixXd priCov = fullCov.block(0, 0, quad.primaryDims, quad.primaryDims);
	Eigen::VectorXd sVar;
	pquad.setup(grp.qwidth, grp.qpoints, grp.mean, priCov, sVar);

	SEXP Rout;
	Rf_protect(Rout = Rf_allocMatrix(REALSXP, outRows, outCols));
	Rf_setAttrib(Rout, R_DimNamesSymbol, dimnames);
	double *out = REAL(Rout);
	memcpy(out, ssProb.data(), sizeof(double) * outRows);
	for (int cx=0; cx <= curMax; cx++) {
		std::vector<double> pad(quad.primaryDims + triangleLoc1(quad.primaryDims));
		pquad.EAP(&slCur.coeffRef(0, cx), 1/ssProb[cx], pad.data());
		for (int sx=0; sx < quad.primaryDims; ++sx) {
			out[(1+sx) * outRows + cx] = pad[sx];
			out[(1+quad.primaryDims+sx) * outRows + cx] = sqrt(pad[quad.primaryDims + triangleLoc0(sx)]);
		}
		for (int sx=0; sx < triangleLoc1(quad.primaryDims); ++sx) {
			out[(1+2*quad.primaryDims+sx) * outRows + cx] = pad[quad.primaryDims + sx];
		}
	}

	return Rout;
}

SEXP pairwiseExpected(SEXP robj, SEXP Rwidth, SEXP Rpts, SEXP Ritems, SEXP Rtwotier)
{
	ProtectAutoBalanceDoodad mpi;

	if (Rf_length(Ritems) != 2) Rf_error("A pair of items must be specified");

	ifaGroup grp(1, Rf_asLogical(Rtwotier));
	grp.setGridFineness(Rf_asReal(Rwidth), Rf_asInteger(Rpts));
	grp.import(robj, false); // lenient=true is probably okay, need to test
	grp.setupQuadrature();
	
	ba81NormalQuad &quad = grp.quad;

	int i1 = INTEGER(Ritems)[0];
	int i2 = INTEGER(Ritems)[1];
	if (i1 < 0 || i1 >= (int) grp.spec.size()) Rf_error("Item %d out of range", i1);
	if (i2 < 0 || i2 >= (int) grp.spec.size()) Rf_error("Item %d out of range", i2);
	if (i1 == i2) Rf_warning("Request to create bivariate distribution of %d with itself", i1);

	double *i1par = &grp.param[i1 * grp.paramRows];
	double *i2par = &grp.param[i2 * grp.paramRows];

	int specific1 = -1;
	int specific2 = -1;
	if (grp.numSpecific) {
		int priDims = quad.maxDims-1;
		for (int ax=priDims; ax < quad.maxAbilities; ax++) {
			if (i1par[ax] != 0) {
				specific1 = ax - priDims;
			}
			if (i2par[ax] != 0) {
				specific2 = ax - priDims;
			}
		}
	}

	const double *spec1 = grp.spec[i1];
	int id1 = spec1[RPF_ISpecID];
	int outcomes1 = spec1[RPF_ISpecOutcomes];

	const double *spec2 = grp.spec[i2];
	int id2 = spec2[RPF_ISpecID];
	int outcomes2 = spec2[RPF_ISpecOutcomes];

	SEXP Rexpected;
	Rf_protect(Rexpected = Rf_allocMatrix(REALSXP, outcomes1, outcomes2));
	double *outMem = REAL(Rexpected);
	Eigen::Map<Eigen::MatrixXd> out(outMem, outcomes1, outcomes2);
	out.setZero();

	// See Cai & Hansen (2012) Eqn 25, 26

	Eigen::VectorXd o1(outcomes1);
	Eigen::VectorXd o2(outcomes2);

	if (specific1 == -1 && specific2 == -1) {
		int specificIncr = quad.numSpecific? quad.quadGridSize : 1;
		for (int qx=0; qx < quad.totalPrimaryPoints; ++qx) {
			double *where = quad.wherePrep.data() + qx * quad.maxDims * specificIncr;
			(*Glibrpf_model[id1].prob)(spec1, i1par, where, o1.data());
			(*Glibrpf_model[id2].prob)(spec2, i2par, where, o2.data());
			out += (o1 * o2.transpose()) * quad.priQarea[qx];
		}
	} else if (specific1 == specific2) {
		Eigen::VectorXd ptheta(quad.maxAbilities);
		for (int qloc=0, qx=0; qx < quad.totalPrimaryPoints; ++qx) {
			for (int sx=0; sx < quad.quadGridSize; ++sx) {
				double *where = quad.wherePrep.data() + qloc * quad.maxDims;
				for (int dx=0; dx < quad.maxAbilities; dx++) {
					ptheta[dx] = where[std::min(dx, quad.maxDims-1)];
				}
				(*Glibrpf_model[id1].prob)(spec1, i1par, ptheta.data(), o1.data());
				(*Glibrpf_model[id2].prob)(spec2, i2par, ptheta.data(), o2.data());
				double area = quad.priQarea[qx] * quad.speQarea[sx * quad.numSpecific + specific1];
				out += (o1 * o2.transpose()) * area;
				++qloc;
			}
		}
	} else if (specific1 != specific2) {
		Eigen::VectorXd spo1(outcomes1);
		Eigen::VectorXd spo2(outcomes2);
		Eigen::VectorXd ptheta(quad.maxAbilities);
		for (int qloc=0, qx=0; qx < quad.totalPrimaryPoints; ++qx) {
			o1.setZero();
			o2.setZero();
			for (int sx=0; sx < quad.quadGridSize; ++sx) {
				double *where = quad.wherePrep.data() + qloc * quad.maxDims;
				for (int dx=0; dx < quad.maxAbilities; dx++) {
					ptheta[dx] = where[std::min(dx, quad.maxDims-1)];
				}
				(*Glibrpf_model[id1].prob)(spec1, i1par, ptheta.data(), spo1.data());
				(*Glibrpf_model[id2].prob)(spec2, i2par, ptheta.data(), spo2.data());
				if (specific1 == -1) {
					if (sx==0) o1 = spo1;
				} else {
					o1 += spo1 * quad.speQarea[sx * quad.numSpecific + specific1];
				}
				if (specific2 == -1) {
					if (sx==0) o2 = spo2;
				} else {
					o2 += spo2 * quad.speQarea[sx * quad.numSpecific + specific2];
				}
				++qloc;
			}
			out += (o1 * o2.transpose()) * quad.priQarea[qx];
		}
	}

	return Rexpected;
}

static const double KANG_CHEN_MIN_EXPECTED = 1.0;  // customizable parameter?

struct ManhattenCollapse {
	Eigen::Map<Eigen::ArrayXXd> obs;
	Eigen::Map<Eigen::ArrayXXd> expected;

	Eigen::DenseIndex smr, smc;
	double bestFit;
	Eigen::DenseIndex bestR, bestC;
	
	ManhattenCollapse(int rows, int cols, double *oMem, double *eMem)
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

SEXP collapse_wrapper(SEXP r_observed_orig, SEXP r_expected_orig)
{
	ProtectAutoBalanceDoodad mpi;

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
  Rf_protect(r_observed = Rf_duplicate(r_observed_orig));
  Rf_protect(r_expected = Rf_duplicate(r_expected_orig));

  double *observed = REAL(r_observed);
  double *expected = REAL(r_expected);

  ManhattenCollapse mcollapse(rows, cols, observed, expected);
  int collapsed = mcollapse.run();

  SEXP dimnames;
  Rf_protect(dimnames = Rf_getAttrib(r_observed_orig, R_DimNamesSymbol));
  Rf_setAttrib(r_observed, R_DimNamesSymbol, dimnames);
  Rf_setAttrib(r_expected, R_DimNamesSymbol, dimnames);

  MxRList out;
  out.add("O", r_observed);
  out.add("E", r_expected);
  out.add("collapsed", Rf_ScalarInteger(collapsed));
  return out.asR();
}

static int maxObservedSumScore(ifaGroup &grp, int *itemMask)
{
	int curMax = 0;
	for (int ix=0; ix < int(grp.spec.size()); ++ix) {
		if (!itemMask[ix]) continue;
		const double *spec = grp.spec[ix];
		int no = spec[RPF_ISpecOutcomes];
		curMax += no - 1;
	}
	return curMax;
}

static bool computeObservedSumScore(ifaGroup &grp, int *itemMask, int row, int *sumOut)
{
	int sum = 0;
	for (int ix=0; ix < int(grp.spec.size()); ++ix) {
		if (!itemMask[ix]) continue;
		const int *resp = grp.dataColumn(ix);
		if (resp[row] == NA_INTEGER) return true;
		sum += resp[row] - 1;
	}
	*sumOut = sum;
	return false;
}

SEXP fast_tableWithWeights(SEXP Ritem1, SEXP Ritem2, SEXP Rweight)
{
	ProtectAutoBalanceDoodad mpi;

	int rows = Rf_length(Ritem1);
	if (rows != Rf_length(Ritem2)) Rf_error("Data are of different lengths");

	Eigen::Map<Eigen::ArrayXi> item1(INTEGER(Ritem1), rows);
	Eigen::Map<Eigen::ArrayXi> item2(INTEGER(Ritem2), rows);
	double *wvec = 0;
	if (Rf_length(Rweight) == rows) wvec = REAL(Rweight);

	SEXP lev1, lev2;
	Rf_protect(lev1 = Rf_getAttrib(Ritem1, R_LevelsSymbol));
	Rf_protect(lev2 = Rf_getAttrib(Ritem2, R_LevelsSymbol));
	int nlev1 = Rf_length(lev1);
	int nlev2 = Rf_length(lev2);

	SEXP Rdist;
	Rf_protect(Rdist = Rf_allocMatrix(REALSXP, nlev1, nlev2));
	Eigen::Map<Eigen::ArrayXXd> result(REAL(Rdist), nlev1, nlev2);
	result.setZero();

	for (int rx=0; rx < rows; ++rx) {
		if (item1[rx] == NA_INTEGER || item2[rx] == NA_INTEGER) continue;
		int i1 = item1[rx] - 1;
		int i2 = item2[rx] - 1;
		double weight = 1;
		if (wvec) weight = wvec[rx];
		result(i1,i2) += weight;
	}

	return Rdist;
}

SEXP observedSumScore(SEXP Rgrp, SEXP Rmask)
{
	ProtectAutoBalanceDoodad mpi;

	ifaGroup grp(1, false);
	grp.import(Rgrp, true);
	if (grp.getNumUnique() == 0) Rf_error("observedSumScore requires data");

	if (Rf_length(Rmask) != int(grp.spec.size())) {
		Rf_error("Mask must be of length %d not %d", int(grp.spec.size()), Rf_length(Rmask));
	}
	int *itemMask = LOGICAL(Rmask);

	int numScores = 1+maxObservedSumScore(grp, itemMask);

	SEXP Rdist;
	Rf_protect(Rdist = Rf_allocVector(REALSXP, numScores));
	Eigen::Map<Eigen::ArrayXd> distOut(REAL(Rdist), numScores);
	distOut.setZero();

	double rowsIncluded = 0;
	for (int rx=0; rx < grp.getNumUnique(); ++rx) {
		int ss;
		if (computeObservedSumScore(grp, itemMask, rx, &ss)) continue;
		double weight = grp.rowWeight? grp.rowWeight[rx] : 1;
		distOut[ss] += weight;
		rowsIncluded += weight;
	}

	MxRList out;
	out.add("dist", Rdist);
	out.add("n", Rf_ScalarReal(rowsIncluded));
	return out.asR();
}

SEXP itemOutcomeBySumScore(SEXP Rgrp, SEXP Rmask, SEXP Rinterest)
{
	ProtectAutoBalanceDoodad mpi;

	ifaGroup grp(1, false);
	grp.import(Rgrp, true);
	if (grp.getNumUnique() == 0) Rf_error("itemOutcomeBySumScore requires data");

	if (Rf_length(Rmask) != int(grp.spec.size())) {
		Rf_error("Mask must be of length %d not %d", int(grp.spec.size()), Rf_length(Rmask));
	}
	int *itemMask = LOGICAL(Rmask);

	int numScores = 1+maxObservedSumScore(grp, itemMask);

	int interest = Rf_asInteger(Rinterest) - 1;
	if (interest < 0 || interest >= int(grp.spec.size())) {
		Rf_error("Item of interest %d must be between 1 and %d", 1+interest, int(grp.spec.size()));
	}

	const double *spec = grp.spec[interest];
	int outcomes = spec[RPF_ISpecOutcomes];

	SEXP r_ans;
	Rf_protect(r_ans = Rf_allocMatrix(REALSXP, numScores, outcomes));
	Eigen::Map<Eigen::ArrayXXd> out(REAL(r_ans), numScores, outcomes);
	out.setZero();

	const int *iresp = grp.dataColumn(interest);

	double rowsIncluded = 0;
	for (int rx=0; rx < grp.getNumUnique(); ++rx) {
		int pick = iresp[rx];
		if (pick == NA_INTEGER) continue;
		int ss;
		if (computeObservedSumScore(grp, itemMask, rx, &ss)) continue;
		double weight = grp.rowWeight? grp.rowWeight[rx] : 1;
		out(ss, pick-1) += weight;
		rowsIncluded += weight;
	}

	MxRList lout;
	lout.add("table", r_ans);
	lout.add("n", Rf_ScalarReal(rowsIncluded));
	return lout.asR();
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
				int r1 = RAND_MAX * unif_rand();
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

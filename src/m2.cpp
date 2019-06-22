#include "rpf.h"

// M2 is not implemented yet

struct ch2012 {
	ifaGroup grp;
	bool pearson;
	double stat;
	double weightSum;
	std::vector<bool> rowMask;

	ch2012(bool twotier, SEXP Rgrp);
	void run(const char *method);
	void accumulate(double observed, double expected);
};

ch2012::ch2012(bool twotier, SEXP Rgrp)
	: grp(1, twotier)
{
	grp.import(Rgrp, false);
	grp.setupQuadrature();
	rowMask.reserve(grp.getNumUnique());
	for (int rx=0; rx < grp.getNumUnique(); ++rx) {
		bool missing = false;
		for (int cx=0; cx < (int) grp.dataColumns.size(); ++cx) {
			if (grp.dataColumns[cx][rx] == NA_INTEGER) {
				missing = true;
				break;
			}
		}
		rowMask.push_back(!missing);
	}
}

void ch2012::accumulate(double observed, double expected)
{
	if (pearson) {
		double diff = observed-expected;
		stat += (diff*diff) / expected;
	} else {
		stat += 2 * observed * (log(observed) - log(expected));
	}
	R_CheckUserInterrupt();  // could loop for a long time
}

void ch2012::run(const char *method)
{
	/*
	std::vector<int> &itemOutcomes = grp.itemOutcomes;

	int numFirstOrder = 0;
	for (int ix=0; ix < grp.numItems(); ++ix) {
		numFirstOrder += itemOutcomes[ix] - 1;
	}

	int numSecondOrder = 0;
	for (int i1=1; i1 < grp.numItems(); ++i1) {
		for (int i2=0; i2 < i1; ++i2) {
			numSecondOrder += (itemOutcomes[i1] - 1) * (itemOutcomes[i2] - 1);
		}
	}
	*/

	if (strEQ(method, "pearson")) {
		pearson = true;
	} else if (strEQ(method, "lr")) {
		pearson = false;
	} else {
		Rf_error("Unknown method '%s'", method);
	}

	if (!grp.rowWeight) Rf_error("weightColumn required");
	ba81NormalQuad &quad = grp.quad;

	grp.ba81OutcomeProb(grp.param, false);

	weightSum = 0;
	for (int rx=0; rx < grp.getNumUnique(); ++rx) {
		if (!rowMask[rx]) continue;
		weightSum += grp.rowWeight[rx];
	}

	stat = 0;
	if (grp.numSpecific == 0) {
		Eigen::ArrayXd Qweight(quad.totalQuadPoints);
		for (int px=0; px < grp.getNumUnique(); ++px) {
			if (!rowMask[px]) continue;
			grp.ba81LikelihoodSlow2(px, Qweight.data());
			accumulate(grp.rowWeight[px], Qweight.sum() * weightSum);
		}
	} else {
		Eigen::ArrayXd Qweight(quad.totalQuadPoints * grp.numSpecific);
		Eigen::ArrayXd Ei(quad.totalPrimaryPoints);
		Eigen::ArrayXd Eis(quad.totalPrimaryPoints * grp.numSpecific);
		for (int px=0; px < grp.getNumUnique(); ++px) {
			if (!rowMask[px]) continue;
			grp.cai2010EiEis(px, Qweight.data(), Eis.data(), Ei.data());
			accumulate(grp.rowWeight[px], Ei.sum() * weightSum);
		}
	}
}

SEXP CaiHansen2012(SEXP Rgrp, SEXP Rmethod, SEXP Rtwotier)
{
	ProtectAutoBalanceDoodad mpi;

	ch2012 engine(Rf_asLogical(Rtwotier), Rgrp);
	engine.run(R_CHAR(Rf_asChar(Rmethod)));
	
	//obMargin1(col);
	//obMargin2(col1, col2);

	MxRList out;
	out.add("stat", Rf_ScalarReal(engine.stat));
	out.add("n", Rf_ScalarReal(engine.weightSum));
	return out.asR();
}

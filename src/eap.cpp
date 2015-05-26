#include "rpf.h"

struct eap {
	std::vector<double *> scoresOut;
};

template <typename T>
struct BA81LatentScores {
       int numLatents;
       Eigen::ArrayXXd thrScore;

       void begin(class ifaGroup *state, T extraData);
       void normalizeWeights(class ifaGroup *state, T extraData, int px, double *Qweight, double weight, int thrId);
       void end(class ifaGroup *state, T extraData);
       bool hasEnd() { return true; }
 };

template <typename T>
void BA81LatentScores<T>::begin(class ifaGroup *state, T extraData)
{
	ba81NormalQuad &quad = state->quad;
	numLatents = quad.maxAbilities + triangleLoc1(quad.maxAbilities);
	thrScore.resize(numLatents, GlobalNumberOfCores);
}

template <typename T>
void BA81LatentScores<T>::normalizeWeights(class ifaGroup *state, T extraData,
					   int px, double *Qweight, double patternLik1, int thrId)
{
	ba81NormalQuad &quad = state->quad;
	const int maxAbilities = quad.maxAbilities;

	// NOTE: Qweight remains unnormalized

	thrScore.col(thrId).setZero();
	double *scorePad = &thrScore.coeffRef(0, thrId);

	quad.EAP(Qweight, 1/patternLik1, scorePad);

	std::vector<double*> &out = extraData.scoresOut;

	for (int ax=0; ax < maxAbilities; ++ax) {
		out[ax][px] = scorePad[ax];
	}
	for (int ax=0; ax < maxAbilities; ++ax) {
		out[maxAbilities + ax][px] = sqrt(scorePad[maxAbilities + triangleLoc0(ax)]);
	}
	for (int ax=0; ax < triangleLoc1(maxAbilities); ++ax) {
		out[2*maxAbilities + ax][px] = scorePad[maxAbilities + ax];
	}
}

template <typename T>
void BA81LatentScores<T>::end(class ifaGroup *state, T extraData)
{
	std::vector<int> &rowMap = state->rowMap;
	const int numUnique = (int) rowMap.size();
	std::vector<double*> &out = extraData.scoresOut;

	for (int px=0; px < numUnique; px++) {
		if (state->patternLik[px]) continue;
		for (int ax=0; ax < int(out.size()); ++ax) {
			out[ax][px] = NA_REAL;
		}
	}
}

SEXP eap_wrapper(SEXP Rgrp)
{
	omxManageProtectInsanity mpi;

	eap eapContext;

	ifaGroup grp(GlobalNumberOfCores, true);
	grp.import(Rgrp);
	grp.buildRowSkip();
	if (grp.getNumUnique() == 0) {
		Rf_error("EAP requested but there are no data rows");
	}
	grp.ba81OutcomeProb(grp.param, false);

	// TODO Wainer & Thissen. (1987). Estimating ability with the wrong
	// model. Journal of Educational Statistics, 12, 339-368.

	/*
	int numQpoints = state->targetQpoints * 2;  // make configurable TODO

	if (numQpoints < 1 + 2.0 * sqrt(state->itemSpec->cols)) {
		// Thissen & Orlando (2001, p. 136)
		Rf_warning("EAP requires at least 2*sqrt(items) quadrature points");
	}

	ba81SetupQuadrature(oo, numQpoints, 0);
	ba81Estep1(oo);
	*/

	int maxAbilities = grp.quad.maxAbilities;
	if (maxAbilities == 0) Rf_error("At least 1 factor is required");
	int rows = grp.getNumUnique();  // allow indexvector for compressed tables TODO
	int cols = 2 * maxAbilities + triangleLoc1(maxAbilities);

	SEXP Rscores;
	Rf_protect(Rscores = Rf_allocVector(VECSXP, cols));
	for (int cx=0; cx < cols; ++cx) {
		SEXP vec = Rf_allocVector(REALSXP, rows);
		SET_VECTOR_ELT(Rscores, cx, vec);
		eapContext.scoresOut.push_back(REAL(vec));
	}

	const int SMALLBUF = 10;
	char buf[SMALLBUF];
	SEXP names;
	Rf_protect(names = Rf_allocVector(STRSXP, cols));
	for (int nx=0; nx < maxAbilities; ++nx) {
		SET_STRING_ELT(names, nx, Rf_mkChar(grp.factorNames[nx].c_str()));
		snprintf(buf, SMALLBUF, "se%d", nx+1);
		SET_STRING_ELT(names, maxAbilities + nx, Rf_mkChar(buf));
	}
	for (int nx=0; nx < triangleLoc1(maxAbilities); ++nx) {
		snprintf(buf, SMALLBUF, "cov%d", nx+1);
		SET_STRING_ELT(names, maxAbilities*2 + nx, Rf_mkChar(buf));
	}
	Rf_setAttrib(Rscores, R_NamesSymbol, names);

	SEXP classes;
	Rf_protect(classes = Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(classes, 0, Rf_mkChar("data.frame"));
	Rf_setAttrib(Rscores, R_ClassSymbol, classes);

	if (grp.dataRowNames) {
		Rf_setAttrib(Rscores, R_RowNamesSymbol, grp.dataRowNames);
	}

	if (grp.numSpecific == 0) {
		BA81Engine<eap&, BA81Dense, BA81LatentScores, BA81OmitEstep> engine;
		engine.ba81Estep1(&grp, eapContext);
	} else {
		BA81Engine<eap&, BA81TwoTier, BA81LatentScores, BA81OmitEstep> engine;
		engine.ba81Estep1(&grp, eapContext);
	}

	return Rscores;
}

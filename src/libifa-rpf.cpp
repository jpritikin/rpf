/*
  Copyright 2012-2014 Joshua Nathaniel Pritikin and contributors

  libifa-rpf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <R.h>
#include <R_ext/BLAS.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../inst/include/libifa-rpf.h"
#include "Eigen/Core"

#ifndef M_LN2
#define M_LN2           0.693147180559945309417232121458        /* ln(2) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI    0.572364942924700087071713675677        /* log(sqrt(pi))
                                                                   == log(pi)/2 */
#endif

extern const struct rpf librpf_model[];
extern const int librpf_numModels;

static const double EXP_STABLE_DOMAIN = 35;
static const double SMALLEST_PROB = 6.305116760146989222002e-16;  // exp(-35), need constexpr

static void
irt_rpf_logprob_adapter(const double *spec,
			const double *param, const double *th,
			double *out)
{
  (*librpf_model[(int) spec[RPF_ISpecID]].prob)(spec, param, th, out);

  int numOutcomes = spec[RPF_ISpecOutcomes];
  for (int ox=0; ox < numOutcomes; ox++) {
    out[ox] = log(out[ox]);
  }
}

static double
dotprod(const double *v1, const double *v2, const int len)
{
  double dprod = 0;
  for (int dx=0; dx < len; dx++) {
    dprod += v1[dx] * v2[dx];
  }
  return dprod;
}

static int
hessianIndex(int numParam, int row, int col)
{
  return numParam + row*(row+1)/2 + col;
}

static double antilogit(const double x)
{
    if (x == INFINITY) return 1.0;
    else if(x == -INFINITY) return 0.0;
    else return 1.0 / (1.0 + exp(-x));
}

static int
irt_rpf_1dim_drm_numSpec(const double *spec)
{ return RPF_ISpecCount; }

static int
irt_rpf_1dim_drm_numParam(const double *spec)
{ return 4; }

static void
irt_rpf_1dim_drm_prob(const double *spec,
		      const double *param, const double *th,
		      double *out)
{
  double guessing = param[2];
  double upper = param[3];
  double athb = -param[0] * (th[0] - param[1]);
  if (athb < -EXP_STABLE_DOMAIN) athb = -EXP_STABLE_DOMAIN;
  else if (athb > EXP_STABLE_DOMAIN) athb = EXP_STABLE_DOMAIN;
  double pp = guessing + (upper-guessing) / (1 + exp(athb));
  out[0] = 1-pp;
  out[1] = pp;
}

static void
set_deriv_nan(const double *spec, double *out)
{
  int numParam = (*librpf_model[(int) spec[RPF_ISpecID]].numParam)(spec);

  for (int px=0; px < numParam; px++) {
    out[px] = nan("I");
  }
}

static void
irt_rpf_1dim_drm_rescale(const double *spec, double *param, const int *paramMask,
			 const double *mean, const double *cov)
{
  double thresh = param[1] * -param[0];
  if (paramMask[0] >= 0) {
    param[0] *= cov[0];
  }
  if (paramMask[1] >= 0) {
    thresh += param[0] * mean[0];
    param[1] = thresh / -param[0];
  }
}

static int
irt_rpf_mdim_drm_numSpec(const double *spec)
{ return RPF_ISpecCount; }

static int
irt_rpf_mdim_drm_numParam(const double *spec)
{
	if (spec[RPF_ISpecDims] == 0) return 1;
	else return 3 + spec[RPF_ISpecDims];
}

static void
irt_rpf_mdim_drm_paramInfo(const double *spec, const int param,
			   const char **type, double *upper, double *lower)
{
	int numDims = spec[RPF_ISpecDims];
	*upper = nan("unset");
	*lower = nan("unset");
	if (numDims == 0) {
		*type = "intercept";
		return;
	}
	*type = NULL;
	if (param >= 0 && param < numDims) {
		*type = "slope";
		*lower = 1e-6;
	} else if (param == numDims) {
		*type = "intercept";
	} else if (param == numDims+1 || param == numDims+2) {
		*type = "bound";
	}
}

static void
irt_rpf_mdim_drm_prob(const double *spec,
		      const double *param, const double *th,
		      double *out)
{
  int numDims = spec[RPF_ISpecDims];
  double dprod = dotprod(param, th, numDims);
  double diff = param[numDims];
  double athb = -(dprod + diff);
  if (athb < -EXP_STABLE_DOMAIN) athb = -EXP_STABLE_DOMAIN;
  else if (athb > EXP_STABLE_DOMAIN) athb = EXP_STABLE_DOMAIN;
  double tmp;
  if (numDims == 0) {
	  tmp = 1/(1+exp(athb));
  } else {
	  const double gg = antilogit(param[numDims+1]);
	  const double uu = antilogit(param[numDims+2]);
	  const double width = uu-gg;
	  if (width < 0) tmp = nan("I");
	  else {
		  tmp = gg + width / (1 + exp(athb));
	  }
  }
  out[0] = 1-tmp;
  out[1] = tmp;
}

static void
irt_rpf_mdim_drm_prob2(const double *spec,
		       const double *param, const double *th,
		       double *out1, double *out2)
{
  int numDims = spec[RPF_ISpecDims];
  double dprod = dotprod(param, th, numDims);
  double diff = param[numDims];
  double athb = -(dprod + diff);
  if (athb < -EXP_STABLE_DOMAIN) athb = -EXP_STABLE_DOMAIN;
  else if (athb > EXP_STABLE_DOMAIN) athb = EXP_STABLE_DOMAIN;
  double tmp = 1 / (1 + exp(athb));
  out1[0] = 1-tmp;
  out1[1] = tmp;
  if (numDims) {
	  const double gg = antilogit(param[numDims+1]);
	  const double uu = antilogit(param[numDims+2]);
	  tmp = gg + (uu-gg) * tmp;
  }
  out2[0] = 1-tmp;
  out2[1] = tmp;
}

static void
irt_rpf_mdim_drm_deriv1(const double *spec,
		       const double *param,
		       const double *where,
		       const double *weight, double *out)
{
  const int numDims = spec[RPF_ISpecDims];
  double QP[2];
  double QPstar[2];
  irt_rpf_mdim_drm_prob2(spec, param, where, QPstar, QP);
  const double r1 = weight[1];
  const double r2 = weight[0];
  const double r1_P = r1/QP[1];
  const double r1_P2 = r1/(QP[1] * QP[1]);
  const double r2_Q = r2/QP[0];
  const double r2_Q2 = r2/(QP[0] * QP[0]);
  const double r1_Pr2_Q = (r1_P - r2_Q);
  const double Pstar = QPstar[1];
  const double Pstar2 = Pstar * Pstar;
  const double Pstar3 = Pstar2 * Pstar;
  const double Qstar = QPstar[0];
  if (numDims == 0) {
	  out[0] -= Pstar * Qstar * r1_Pr2_Q;
	  double chunk1 = (Pstar - Pstar2);
	  out[1] -= (r1_P * ((Pstar - 3*Pstar2 + 2*Pstar3)) -
		     r1_P2 * chunk1*chunk1 +
		     r2_Q * ((-Pstar + 3*Pstar2 - 2*Pstar3)) -
		     r2_Q2 * chunk1*chunk1);   // cc^2
	  return;
  }
  const double expgg = param[numDims+1];
  const double expuu = param[numDims+2];
  const double gg = antilogit(expgg);
  const double uu = antilogit(expuu);
  const double difexpgg = gg * (1-gg);
  const double difexpuu = uu * (1-uu);
  const double gm1 = (1.0 - gg);
  const double um1 = (1.0 - uu);
  const double u_1u = uu * um1;
  const double g_1g = gg * gm1;
  const double ugD = (uu-gg);
  for (int dx=0; dx < numDims; dx++) {
    out[dx] -= where[dx] * Pstar * Qstar * ugD * r1_Pr2_Q;
  }
  out[numDims] -= ugD * Pstar * Qstar * r1_Pr2_Q;
  out[numDims+1] -= difexpgg * QPstar[0] * r1_Pr2_Q;
  out[numDims+2] -= difexpuu * QPstar[1] * r1_Pr2_Q;

  int ox = numDims+2;

  for(int ix=0; ix < numDims; ix++) {
    for(int jx=0; jx <= ix; jx++) {
      out[++ox] -= (r1_P * (ugD * where[ix] * where[jx] *
				   (Pstar - 3*Pstar2 + 2*Pstar3)) -
			   r1_P2 * (ugD * where[ix] * (Pstar - Pstar2) *
				    (ugD * where[jx] * (Pstar - Pstar2))) +
			   r2_Q * (ugD * where[ix] * where[jx] *
				   (-Pstar + 3*Pstar2 - 2*Pstar3)) -
			   r2_Q2 * (ugD * where[ix] * (-Pstar + Pstar2) *
				    (ugD * where[jx] * (-Pstar + Pstar2))));  // aa_k aa_k
    }
  }
  for(int ix=0; ix < numDims; ix++) {
    out[++ox] -= (r1_P * (ugD * where[ix] * (Pstar - 3*Pstar2 + 2*Pstar3)) -
			 r1_P2 * (ugD * where[ix] * (Pstar - Pstar2) *
				  (ugD * (Pstar - Pstar2))) +
			 r2_Q * (ugD * where[ix] * (-Pstar + 3*Pstar2 - 2*Pstar3)) -
			 r2_Q2 * (ugD * where[ix] * (-Pstar + Pstar2) *
				  (ugD * (-Pstar + Pstar2))));  // cc aa_k
  }
  double chunk1 = ugD * (Pstar - Pstar2);
  out[++ox] -= (r1_P * (ugD * (Pstar - 3*Pstar2 + 2*Pstar3)) -
		r1_P2 * chunk1*chunk1 +
		r2_Q * (ugD * (-Pstar + 3*Pstar2 - 2*Pstar3)) -
		r2_Q2 * chunk1*chunk1);   // cc^2
  for(int ix=0; ix < numDims; ix++) {
	  out[++ox] -= (r1_P * (g_1g * where[ix] * (-Pstar + Pstar2)) -
		  r1_P2 * (ugD * where[ix] * (Pstar - Pstar2)) * g_1g * Qstar +
		  r2_Q * (g_1g * where[ix] * (Pstar - Pstar2)) -
		  r2_Q2 * (ugD * where[ix] * (-Pstar + Pstar2) ) * g_1g * (Pstar - 1));   // gg aa_k
  }
  out[++ox] -= (r1_P * (g_1g * (-Pstar + Pstar2)) -
		r1_P2 * (ugD * (Pstar - Pstar2)) * g_1g * Qstar +
		r2_Q * (g_1g * (Pstar - Pstar2)) -
		r2_Q2 * (ugD * (-Pstar + Pstar2)) * g_1g * -Qstar);  // gg cc
  out[++ox] -= (r1_P * (g_1g * (2.0*gm1 - 1.0 - 2.0*gm1*Pstar + Pstar)) -
		r1_P2 * (g_1g * (1.0 - Pstar)) * (g_1g * (1.0 - Pstar)) +
		r2_Q * (g_1g * (-2.0*gm1 + 1.0 + 2.0*gm1*Pstar - Pstar)) -
		r2_Q2 * (g_1g * (-1.0 + Pstar)) * (g_1g * (-1.0 + Pstar)));  // gg^2

  for(int ix=0; ix < numDims; ix++) {
    out[++ox] -= (r1_P * (u_1u * where[ix] * (Pstar - Pstar2)) -
		  r1_P2 * (ugD * where[ix] * (Pstar - Pstar2)) * u_1u * Pstar +
		  r2_Q * (u_1u * where[ix] * (-Pstar + Pstar2)) +
		  r2_Q2 * (ugD * where[ix] * (-Pstar + Pstar2) ) * u_1u * Pstar);  // uu aa_k
  }

  out[++ox] -= (r1_P * (u_1u * (Pstar - Pstar2)) -
		r1_P2 * (ugD * (Pstar - Pstar2)) * u_1u * Pstar +
		r2_Q * (u_1u * (-Pstar + Pstar2)) +
		r2_Q2 * (ugD * (-Pstar + Pstar2)) * u_1u * Pstar);  // uu cc

  out[++ox] -= (-r1_P2 * (g_1g * (1.0 - Pstar)) * u_1u * Pstar +
		r2_Q2 * (g_1g * (-1.0 + Pstar)) * u_1u * Pstar);  // uu gg
  out[++ox] -=  (r1_P * (2.0*u_1u*um1*Pstar) - r1_P * (u_1u*Pstar) - r1_P2 *(u_1u*u_1u*Pstar2) -
		 r2_Q * (2.0*u_1u*um1*Pstar) + r2_Q * (u_1u*Pstar) - r2_Q2 *(u_1u*u_1u*Pstar2));  // uu^2
}

static void
irt_rpf_mdim_drm_deriv2(const double *spec,
			const double *param,
			double *out)
{
  int numDims = spec[RPF_ISpecDims];
  if (numDims == 0) return;
  const double *aa = param;
  double gg = param[numDims+1];
  double uu = param[numDims+2];

  for (int dx=0; dx < numDims; dx++) {
    if (aa[dx] < 0) {
      set_deriv_nan(spec, out);
      return;
    }
  }
  if (gg == -INFINITY) {
    out[numDims+1] = nan("I");
  }
  if (uu == INFINITY) {
    out[numDims+2] = nan("I");
  }
  if (gg > uu) {
    out[numDims+1] = nan("I");
    out[numDims+2] = nan("I");
  }
}

static void
irt_rpf_mdim_drm_rescale(const double *spec, double *param, const int *paramMask,
			 const double *mean, const double *cov)
{
  int numDims = spec[RPF_ISpecDims];

  double madj = dotprod(param, mean, numDims);

  for (int d1=0; d1 < numDims; d1++) {
    if (paramMask[d1] < 0) continue;
    param[d1] = dotprod(param+d1, cov + d1 * numDims + d1, numDims-d1);
  }

  param[numDims] += madj;
}

static void
irt_rpf_mdim_drm_dTheta(const double *spec, const double *param,
			const double *where, const double *dir,
			double *grad, double *hess)
{
  int numDims = spec[RPF_ISpecDims];
  double PQ[2];
  double PQstar[2];
  irt_rpf_mdim_drm_prob2(spec, param, where, PQstar, PQ);
  double Pstar = PQstar[0];
  double Qstar = PQstar[1];
  const double *aa = param;
  const double guess = antilogit(param[numDims + 1]);
  const double upper = antilogit(param[numDims + 2]);
  for (int ax=0; ax < numDims; ax++) {
    double piece = dir[ax] * (upper-guess) * aa[ax] * (Pstar * Qstar);
    grad[1] += piece;
    grad[0] -= piece;
    piece = dir[ax] * (2 * (upper - guess) * aa[ax]*aa[ax] * (Qstar * Qstar * Pstar) -
		       (upper - guess) * aa[ax]*aa[ax] * (Pstar * Qstar));
    hess[1] -= piece;
    hess[0] += piece;
  }
}

static void
irt_rpf_1dim_drm_dTheta(const double *spec, const double *param,
			const double *where, const double *dir,
			double *grad, double *hess)
{
  double nparam[4];
  memcpy(nparam, param, sizeof(double) * 4);
  nparam[1] = param[1] * -param[0];
  irt_rpf_mdim_drm_dTheta(spec, nparam, where, dir, grad, hess);
}

static int
irt_rpf_mdim_grm_numSpec(const double *spec)
{ return RPF_ISpecCount; }

static int
irt_rpf_mdim_grm_numParam(const double *spec)
{ return spec[RPF_ISpecOutcomes] + spec[RPF_ISpecDims] - 1; }

static void
irt_rpf_mdim_grm_paramInfo(const double *spec, const int param,
			   const char **type, double *upper, double *lower)
{
	int numDims = spec[RPF_ISpecDims];
	*upper = nan("unset");
	*lower = nan("unset");
	*type = NULL;
	if (param >= 0 && param < numDims) {
		*type = "slope";
		*lower = 1e-6;
	} else {
		*type = "intercept";
	}
}

static void _grm_fix_crazy_stuff(const double *spec, const int numOutcomes, double *out)
{
  int bigk = -1;
  double big = 0;

  for (int bx=0; bx < numOutcomes; bx++) {
    if (out[bx] > big) {
      bigk = bx;
      big = out[bx];
    }
  }

  for (int fx=0; fx < numOutcomes; fx++) {
    if (out[fx] < SMALLEST_PROB) {
      double small = SMALLEST_PROB - out[fx];
      out[bigk] -= small;
      out[fx] += small;
    }
  }
}

static void
irt_rpf_mdim_grm_prob(const double *spec,
		      const double *param, const double *th,
		      double *out)
{
  const int numDims = spec[RPF_ISpecDims];
  const int numOutcomes = spec[RPF_ISpecOutcomes];
  const double *slope = param;
  const double dprod = dotprod(slope, th, numDims);
  const double *kat = param + (int) spec[RPF_ISpecDims];

  double athb = -(dprod + kat[0]);
  if (athb < -EXP_STABLE_DOMAIN) athb = -EXP_STABLE_DOMAIN;
  else if (athb > EXP_STABLE_DOMAIN) athb = EXP_STABLE_DOMAIN;
  double tmp = 1 / (1 + exp(athb));
  out[0] = 1-tmp;
  out[1] = tmp;

  for (int kx=2; kx < numOutcomes; kx++) {
	  if (1e-6 + kat[kx-1] >= kat[kx-2]) {
		  for (int ky=0; ky < numOutcomes; ky++) {
			  out[ky] = nan("I");
		  }
		  return;
	  }
    double athb = -(dprod + kat[kx-1]);
    if (athb < -EXP_STABLE_DOMAIN) athb = -EXP_STABLE_DOMAIN;
    else if (athb > EXP_STABLE_DOMAIN) athb = EXP_STABLE_DOMAIN;
    double tmp = 1 / (1 + exp(athb));
    out[kx-1] -= tmp;
    out[kx] = tmp;
  }

  for (int kx=0; kx < numOutcomes; kx++) {
    if (out[kx] <= 0) {
      _grm_fix_crazy_stuff(spec, numOutcomes, out);
      return;
    }
  }
}

static void
irt_rpf_mdim_grm_rawprob(const double *spec,
			 const double *param, const double *th,
			 double *out)
{
  int numDims = spec[RPF_ISpecDims];
  const int numOutcomes = spec[RPF_ISpecOutcomes];
  const double dprod = dotprod(param, th, numDims);
  const double *kat = param + (int) spec[RPF_ISpecDims];

  out[0] = 1;
  for (int kx=0; kx < numOutcomes-1; kx++) {
    double athb = -(dprod + kat[kx]);
    if (athb < -EXP_STABLE_DOMAIN) athb = -EXP_STABLE_DOMAIN;
    else if (athb > EXP_STABLE_DOMAIN) athb = EXP_STABLE_DOMAIN;
    double tmp = 1 / (1 + exp(athb));
    out[kx+1] = tmp;
  }
  out[numOutcomes] = 0;
}

// Compare with Cai (2010, p. 54) Appendix B
static void
irt_rpf_mdim_grm_deriv1(const double *spec,
			const double *param,
			const double *where,
			const double *weight, double *out)
{
  int nfact = spec[RPF_ISpecDims];
  int outcomes = spec[RPF_ISpecOutcomes];
  int nzeta = spec[RPF_ISpecOutcomes] - 1;
  Eigen::VectorXd P(nzeta+2);
  Eigen::VectorXd PQfull(nzeta+2);
  irt_rpf_mdim_grm_rawprob(spec, param, where, P.data());
  PQfull[0] = 0;
  PQfull[outcomes] = 0;
  for (int kx=1; kx <= nzeta; kx++) PQfull[kx] = P[kx] * (1-P[kx]);
  for (int jx = 0; jx <= nzeta; jx++) {
    double Pk_1 = P[jx];
    double Pk = P[jx + 1];
    double PQ_1 = PQfull[jx];
    double PQ = PQfull[jx + 1];
    double Pk_1Pk = Pk_1 - Pk;
    if (Pk_1Pk < 1e-10) Pk_1Pk = 1e-10;
    double dif1 = weight[jx] / Pk_1Pk;
    double dif1sq = dif1 / Pk_1Pk;
    if(jx < nzeta) {
      double Pk_p1 = P[jx + 2];
      double PQ_p1 = PQfull[jx + 2];
      double Pk_Pkp1 = Pk - Pk_p1;
      if(Pk_Pkp1 < 1e-10) Pk_Pkp1 = 1e-10;
      double dif2 = weight[jx+1] / Pk_Pkp1;
      double dif2sq = dif2 / Pk_Pkp1;
      out[nfact + jx] += PQ * (dif1 - dif2);  //gradient for intercepts

      int d2base = hessianIndex(nfact + nzeta, nfact+jx, 0);
      // hessian for intercept^2
      double tmp3 = (dif1 - dif2) * (Pk * (1.0 - Pk) * (1.0 - 2.0*Pk));
      double piece1 = (PQ * PQ * (dif1sq + dif2sq) + tmp3);
      out[d2base + nfact + jx] += piece1;
      if (jx < (nzeta - 1)) {
	      // hessian for adjacent intercepts
	      int d2base1 = hessianIndex(nfact + nzeta, nfact+jx+1, nfact + jx);
	      out[d2base1] -= dif2sq * PQ_p1 * PQ;
      }
      double tmp1 = -dif2sq * PQ * (PQ - PQ_p1);
      double tmp2 = dif1sq * PQ * (PQ_1 - PQ);
      for(int kx = 0; kx < nfact; kx++){
	// hessian for slope intercept
	out[d2base + kx] -= (tmp1 + tmp2 - tmp3) * where[kx];
      }
    }
    for(int kx = 0; kx < nfact; kx++) {
      // gradient for slope
      out[kx] -= dif1 * (PQ_1 - PQ) * where[kx];
    }

    Eigen::VectorXd temp(nfact);
    for(int ix = 0; ix < nfact; ix++)
      temp[ix] = PQ_1 * where[ix] - PQ * where[ix];

    int d2x = nfact + nzeta;
    double Pk_adj = (Pk_1 * (1.0 - Pk_1) * (1.0 - 2.0 * Pk_1) -
		     Pk * (1.0 - Pk) * (1.0 - 2.0 * Pk));
    for(int i = 0; i < nfact; i++) {
      for(int j = 0; j <= i; j++) {
	double outer = where[i]*where[j];
	// hessian for slope slope
	out[d2x++] -= (- dif1sq * temp[i] * temp[j] + (dif1 * outer * Pk_adj));
      }
    }
  }
}

static void
irt_rpf_mdim_grm_deriv2(const double *spec,
			const double *param,
			double *out)
{
  int nfact = spec[RPF_ISpecDims];
  int nzeta = spec[RPF_ISpecOutcomes] - 1;
  const double *aa = param;
  for (int dx=0; dx < nfact; dx++) {
    if (aa[dx] < 0) {
      set_deriv_nan(spec, out);
      return;
    }
  }
  for (int zx=0; zx < nzeta-1; zx++) {
    if (param[nfact+zx] < param[nfact+zx+1]) {
      set_deriv_nan(spec, out);
      return;
    }
  }
}

static void
irt_rpf_mdim_grm_dTheta(const double *spec, const double *param,
			const double *where, const double *dir,
			double *grad, double *hess)
{
  int numDims = spec[RPF_ISpecDims];
  int outcomes = spec[RPF_ISpecOutcomes];
  const double *aa = param;
  Eigen::VectorXd P(outcomes+1);
  irt_rpf_mdim_grm_rawprob(spec, param, where, P.data());
  for (int jx=0; jx < numDims; jx++) {
    for (int ix=0; ix < outcomes; ix++) {
      double w1 = P[ix] * (1-P[ix]) * aa[jx];
      double w2 = P[ix+1] * (1-P[ix+1]) * aa[jx];
      grad[ix] += dir[jx] * (w1 - w2);
      hess[ix] += dir[jx] * (aa[jx]*aa[jx] * (2 * P[ix] * (1 - P[ix])*(1 - P[ix]) -
					      P[ix] * (1 - P[ix]) -
					      2 * P[ix+1] * (1 - P[ix+1])*(1 - P[ix+1]) +
					      P[ix+1] * (1 - P[ix+1])));
    }
  }
}

static void
irt_rpf_mdim_grm_rescale(const double *spec, double *param, const int *paramMask,
			 const double *mean, const double *cov)
{
  int numDims = spec[RPF_ISpecDims];
  int nzeta = spec[RPF_ISpecOutcomes] - 1;

  double madj = dotprod(param, mean, numDims);

  for (int d1=0; d1 < numDims; d1++) {
    if (paramMask[d1] < 0) continue;
    param[d1] = dotprod(param+d1, cov + d1 * numDims + d1, numDims-d1);
  }

  for (int tx=0; tx < nzeta; tx++) {
    int px = numDims + tx;
    if (paramMask[px] >= 0) param[px] += madj;
  }
}

static int
irt_rpf_nominal_numSpec(const double *spec)
{
  int outcomes = spec[RPF_ISpecOutcomes];
  int Tlen = (outcomes - 1) * (outcomes - 1);
  return RPF_ISpecCount + 4 * Tlen;
}

static int
irt_rpf_nominal_numParam(const double *spec)
{
	int dims = spec[RPF_ISpecDims];
	if (dims == 0) {
		return spec[RPF_ISpecOutcomes]-1;
	} else {
		return dims + 2 * (spec[RPF_ISpecOutcomes]-1);
	}
}


static void
irt_rpf_nominal_paramInfo(const double *spec, const int param,
			  const char **type, double *upper, double *lower)
{
	int numDims = spec[RPF_ISpecDims];
	const int numOutcomes = spec[RPF_ISpecOutcomes];
	*upper = nan("unset");
	*lower = nan("unset");
	if (numDims == 0) {
		*type = "intercept";
		return;
	}
	*type = NULL;
	if (param >= 0 && param < numDims) {
		*type = "slope";
		*lower = 1e-6;
	} else if (param < numDims + numOutcomes - 1) {
		*type = "slope";
	} else {
		*type = "intercept";
	}
}

static void
_nominal_rawprob1(const double *spec,
		 const double *param, const double *th,
		 double discr, double *ak, double *num, double *maxout)
{
  int numDims = spec[RPF_ISpecDims];
  int numOutcomes = spec[RPF_ISpecOutcomes];
  const double *alpha = param + numDims;
  const double *gamma;
  if (numDims == 0) {
	  // alpha doesn't matter
	  gamma = param + numDims;
  } else {
	  gamma = param + numDims + numOutcomes - 1;
  }
  const double *Ta = spec + RPF_ISpecCount;
  const double *Tc = spec + RPF_ISpecCount + (numOutcomes-1) * (numOutcomes-1);

  double curmax = 1;
  for (int kx=0; kx < numOutcomes; kx++) {
    ak[kx] = 0;
    double ck = 0;
    if (kx) {
      for (int tx=0; tx < numOutcomes-1; tx++) {
	int Tcell = tx * (numOutcomes-1) + kx-1;
	ak[kx] += Ta[Tcell] * alpha[tx];
	ck += Tc[Tcell] * gamma[tx];
      }
    }

    double z = discr * ak[kx] + ck;
    num[kx] = z;
    if (curmax < z) curmax = z;
  }
  *maxout = curmax;
}

static void
_nominal_rawprob2(const double *spec,
		  const double *param, const double *th,
		  double discr, double *ak, double *num)
{
  int numOutcomes = spec[RPF_ISpecOutcomes];
  double maxZ;
  _nominal_rawprob1(spec, param, th, discr, ak, num, &maxZ);

  double recenter = 0;
  if (maxZ > EXP_STABLE_DOMAIN) {
    recenter = maxZ - EXP_STABLE_DOMAIN;
  }

  int Kadj = -1;
  double adj = 0;
  double den = 0;   // not exact because adj not taken into account
  for (int kx=0; kx < numOutcomes; kx++) {
    if (num[kx] == maxZ) Kadj = kx;
    if (num[kx] - recenter < -EXP_STABLE_DOMAIN) {
      num[kx] = 0;
      adj += SMALLEST_PROB;
      continue;
    }
    num[kx] = exp(num[kx] - recenter);
    den += num[kx];
  }
  for (int kx=0; kx < numOutcomes; kx++) {
    if (kx == Kadj) {
      num[kx] = num[kx]/den - adj;
    } else if (num[kx] == 0) {
      num[kx] = SMALLEST_PROB;
    } else {
      num[kx] = num[kx]/den;
    }
  }
}

static void
irt_rpf_nominal_prob(const double *spec,
		     const double *param, const double *th,
		     double *out)
{
  int numOutcomes = spec[RPF_ISpecOutcomes];
  int numDims = spec[RPF_ISpecDims];
  Eigen::VectorXd ak(numOutcomes);
  double discr = dotprod(param, th, numDims);
  _nominal_rawprob2(spec, param, th, discr, ak.data(), out);
}

static void
irt_rpf_nominal_logprob(const double *spec,
			const double *param, const double *th,
			double *out)
{
  int numOutcomes = spec[RPF_ISpecOutcomes];
  int numDims = spec[RPF_ISpecDims];
  Eigen::VectorXd num(numOutcomes);
  Eigen::VectorXd ak(numOutcomes);
  double discr = dotprod(param, th, numDims);
  double maxZ;
  _nominal_rawprob1(spec, param, th, discr, ak.data(), num.data(), &maxZ);
  double den = 0;

  if (maxZ > EXP_STABLE_DOMAIN) {
    den = maxZ;  // not best approx
  } else {
    for (int kx=0; kx < numOutcomes; kx++) {
      if (num[kx] < -EXP_STABLE_DOMAIN) continue;
      den += exp(num[kx]);
    }
    den = log(den);
  }

  for (int kx=0; kx < numOutcomes; kx++) {
    out[kx] = num[kx] - den;
  }
}

static double makeOffterm(const double *dat, const double p, const double aTheta,
			  const int ncat, const int cat)
{
  double ret = 0;
  for (int CAT = 0; CAT < ncat; CAT++) {
    if (CAT == cat) continue;
    ret += dat[CAT] * p * aTheta;
  }
  return(ret);
}

static double makeOffterm2(const double *dat, const double p1, const double p2,
			   const double aTheta, const int ncat, const int cat)
{
  double ret = 0;
  for (int CAT = 0; CAT < ncat; CAT++) {
    if (CAT == cat) continue;
    ret += dat[CAT] * p1 * p2 * aTheta;
  }
  return(ret);
}

static void
irt_rpf_nominal_deriv1(const double *spec,
		       const double *param,
		       const double *where,
		       const double *weight, double *out)
{
  int nfact = spec[RPF_ISpecDims];
  int ncat = spec[RPF_ISpecOutcomes];
  double aTheta = dotprod(param, where, nfact);
  double aTheta2 = aTheta * aTheta;

  Eigen::VectorXd num(ncat);
  Eigen::VectorXd ak(ncat);
  _nominal_rawprob2(spec, param, where, aTheta, ak.data(), num.data());

  Eigen::VectorXd P(ncat);
  Eigen::VectorXd P2(ncat);
  Eigen::VectorXd P3(ncat);
  Eigen::VectorXd ak2(ncat);
  Eigen::VectorXd dat_num(ncat);
  double numsum = 0;
  double numakD = 0;
  double numak2D2 = 0;
  Eigen::VectorXd numakDTheta_numsum(nfact);

  for (int kx=0; kx < ncat; kx++) {
    ak2[kx] = ak[kx] * ak[kx];
    dat_num[kx] = weight[kx]/num[kx];
    numsum += num[kx];
    numakD += num[kx] * ak[kx];
    numak2D2 += num[kx] * ak2[kx];
  }
  double numsum2 = numsum * numsum;

  for (int kx=0; kx < ncat; kx++) {
    P[kx] = num[kx]/numsum;
    P2[kx] = P[kx] * P[kx];
    P3[kx] = P2[kx] * P[kx];
  }

  double sumNumak = dotprod(num.data(), ak.data(), ncat);
  for (int fx=0; fx < nfact; fx++) {
    numakDTheta_numsum[fx] = sumNumak * where[fx] / numsum;
  }

  for (int jx = 0; jx < nfact; jx++) {
    double tmpvec = 0;
    for(int i = 0; i < ncat; i++) {
      tmpvec += dat_num[i] * (ak[i] * where[jx] * P[i] -
			      P[i] * numakDTheta_numsum[jx]) * numsum;
    }
    out[jx] -= tmpvec;
  }
  int dkoffset;
  if (nfact == 0) {
	  dkoffset = 0;
  } else {
	  dkoffset = ncat - 1;
  }
  for(int i = 1; i < ncat; i++) {
	  if (nfact) {
		  double offterm = makeOffterm(weight, P[i], aTheta, ncat, i);
		  double tmpvec = dat_num[i] * (aTheta * P[i] - P2[i] * aTheta) * numsum - offterm;
		  out[nfact + i - 1] -= tmpvec;
	  }
    double offterm2 = makeOffterm(weight, P[i], 1, ncat, i);
    double tmpvec2 = dat_num[i] * (P[i] - P2[i]) * numsum - offterm2;
    out[nfact + dkoffset + i - 1] -= tmpvec2;
  }

  int hessbase = nfact + (ncat-1) + dkoffset;
  int d2ind = 0;
  //a's
  for (int j = 0; j < nfact; j++) {
    for (int k = 0; k <= j; k++) {
      double tmpvec = 0;
      for (int i = 0; i < ncat; i++) {
	tmpvec += dat_num[i] * (ak2[i] * where[j] * where[k] * P[i] -
				ak[i] * where[j] * P[i] * numakDTheta_numsum[k] -
				ak[i] * where[k] * P[i] * numakDTheta_numsum[j] +
				2 * P[i] * numakD * where[j] * numakD * where[k] / numsum2 -
				P[i] * numak2D2 * where[j] * where[k] / numsum) * numsum -
	  dat_num[i] * (ak[i] * where[j] * P[i] - P[i] * numakDTheta_numsum[j]) *
	  numsum * ak[i] * where[k] +
	  dat_num[i] * (ak[i] * where[j] * P[i] - P[i] * numakDTheta_numsum[j]) *
	  numakD * where[k];
      }
      out[hessbase + d2ind++] -= tmpvec;
    }
  }
  //a's with ak and d
  for(int k = 1; k < ncat; k++){
    int akrow = hessbase + (nfact+k)*(nfact+k-1)/2;
    int dkrow = hessbase + (nfact+ncat+k-1)*(nfact+ncat+k-2)/2;
    for(int j = 0; j < nfact; j++){
      double tmpvec = 0;
      double tmpvec2 = 0;
      for(int i = 0; i < ncat; i++){
	if(i == k){
	  tmpvec += dat_num[i] * (ak[i]*where[j] * aTheta*P[i] -
				     aTheta*P[i]*numakDTheta_numsum[j] +
				     where[j]*P[i] - 2*ak[i]*where[j]*aTheta*P2[i] +
				     2*aTheta*P2[i]*numakDTheta_numsum[j] -
				     where[j]*P2[i])*numsum -
	    dat_num[i]*(aTheta*P[i] - aTheta*P2[i])*numsum*ak[i]*where[j] +
	    dat_num[i]*(aTheta*P[i] - aTheta*P2[i])*(numakD*where[j]);
	  tmpvec2 += dat_num[i]*(ak[i]*where[j]*P[i] -
				      2*ak[i]*where[j]*P2[i] -
				      P[i]*numakDTheta_numsum[j] +
				      2*P2[i]*numakDTheta_numsum[j])*numsum -
	    dat_num[i]*(P[i] - P2[i])*numsum*ak[i]*where[j] +
	    dat_num[i]*(P[i] - P2[i])*(numakD*where[j]);
	} else {
	  tmpvec += -weight[i]*ak[k]*aTheta*where[j]*P[k] +
	    weight[i]*P[k]*aTheta*numakDTheta_numsum[j] -
	    weight[i]*P[k]*where[j];
	  tmpvec2 += -weight[i]*ak[k]*where[j]*P[k] +
	    weight[i]*P[k]*numakDTheta_numsum[j];
	}
      }
      out[akrow + j] -= tmpvec;
      out[dkrow + j] -= tmpvec2;
    }
  }
  //ak's and d's
  for(int j = 1; j < ncat; j++){
    int akrow = hessbase + (nfact+j)*(nfact+j-1)/2;
    int dkrow = hessbase + (nfact+dkoffset+j)*(nfact+dkoffset+j-1)/2;

    double tmpvec = makeOffterm(weight, P2[j], aTheta2, ncat, j);
    double tmpvec2 = makeOffterm(weight, P[j], aTheta2, ncat, j);
    double offterm = tmpvec - tmpvec2;
    tmpvec = makeOffterm(weight, P2[j], 1, ncat, j);
    tmpvec2 = makeOffterm(weight, P[j], 1, ncat, j);
    double offterm2 = tmpvec - tmpvec2;

    if (nfact) {
	    out[akrow + nfact + j - 1] -=
		    (dat_num[j]*(aTheta2*P[j] - 3*aTheta2*P2[j] +
				 2*aTheta2*P3[j])*numsum - weight[j]/num[j] *
		     (aTheta*P[j] - aTheta*P2[j])*numsum*aTheta + weight[j] *
		     (aTheta*P[j] - aTheta*P2[j])*aTheta + offterm);
    }

    out[dkrow + nfact + dkoffset + j - 1] -=
      (dat_num[j]*(P[j] - 3*P2[j] + 2*P3[j])*numsum - weight[j]/num[j] *
	      (P[j] - P2[j])*numsum + weight[j] *
	      (P[j] - P2[j]) + offterm2);

    for(int i = 1; i < ncat; i++) {
      if(j > i) {
	      if (nfact) {
		      offterm = makeOffterm2(weight, P[j], P[i], aTheta2, ncat, i);
		      tmpvec = dat_num[i] * (-aTheta2*P[i]*P[j] + 2*P2[i] *aTheta2*P[j])*numsum +
			      dat_num[i] * (aTheta*P[i] - P2[i] * aTheta)*aTheta*num[j]+offterm;
		      out[akrow + nfact + i - 1] -= tmpvec;
	      }
	offterm2 = makeOffterm2(weight, P[j], P[i], 1, ncat, i);
	tmpvec2 = dat_num[i] * (-P[i]*P[j] + 2*P2[i] *P[j]) * numsum +
	  dat_num[i] * (P[i] - P2[i]) * num[j] + offterm2;
	out[dkrow + nfact + dkoffset + i - 1] -= tmpvec2;
      }
      if (nfact == 0) continue;
      if (abs(j-i) == 0) {
	tmpvec = makeOffterm(weight, P2[i], aTheta, ncat, i);
	tmpvec2 = makeOffterm(weight, P[i], aTheta, ncat, i);
	offterm = tmpvec - tmpvec2;
	tmpvec = dat_num[i]*(aTheta*P[i] - 3*aTheta*P2[i] +
			     2*aTheta*P3[i]) * numsum - dat_num[i] *
	  (aTheta*P[i] - aTheta*P2[i])*numsum + weight[i] *
	  (P[i] - P2[i])*aTheta + offterm;
	out[dkrow + nfact + i - 1] -= tmpvec;
      } else {
	offterm = makeOffterm2(weight, P[j], P[i], aTheta, ncat, i);
	tmpvec = dat_num[i] * (-aTheta*P[i]*P[j] + 2*P2[i] *aTheta*P[j]) * numsum +
	  dat_num[i] * (P[i] - P2[i]) * aTheta * num[j] + offterm;
	out[dkrow + nfact + i - 1] -= tmpvec;
      }
    }
  }
}

static void
irt_rpf_nominal_deriv2(const double *spec,
		       const double *param,
		       double *out)
{
  int nfact = spec[RPF_ISpecDims];
  int nzeta = spec[RPF_ISpecOutcomes] - 1;
  const double *aa = param;

  for (int dx=0; dx < nfact; dx++) {
    if (aa[dx] < 0) {
      set_deriv_nan(spec, out);
      return;
    }
  }

  int ckoffset = nzeta;
  if (nfact == 0) ckoffset = 0;

  const double *Ta = spec + RPF_ISpecCount;
  const double *Tc = spec + RPF_ISpecCount + nzeta * nzeta;
  const int numParam = irt_rpf_nominal_numParam(spec);
  Eigen::VectorXd rawOut(numParam);
  memcpy(rawOut.data(), out, sizeof(double) * numParam);

  // gradient
  for (int tx=0; tx < nzeta; tx++) {
    double ak1=0;
    double ck1=0;
    for (int kx=0; kx < nzeta; kx++) {
      int Tcell = tx * nzeta + kx;
      ak1 += rawOut[nfact + kx] * Ta[Tcell];
      ck1 += rawOut[nfact + ckoffset + kx] * Tc[Tcell];
    }
    out[nfact + tx] = ak1;
    out[nfact + ckoffset + tx] = ck1;
  }

  // don't need to transform the main a parameters TODO
  double *dmat = Realloc(NULL, 3 * numParam * numParam, double);
  const int hsize = hessianIndex(0, numParam-1, numParam-1);
  {
	  // unpack triangular storage into a full matrix
    int row=0;
    int col=0;
    for (int dx=0; dx <= hsize; dx++) {
      dmat[numParam * col + row] = out[numParam + dx];
      if (row == col) {
	col=0; ++row;
      } else {
	dmat[numParam * row + col] = out[numParam + dx];
	++col;
      }
    }
  }

  double *tmat = dmat + numParam * numParam;
  for (int dx=0; dx < numParam * numParam; dx++) tmat[dx] = 0;
  for (int dx=0; dx < nfact; dx++) {
    tmat[dx * numParam + dx] = 1;
  }
  for (int rx=0; rx < nzeta; rx++) {
    for (int cx=0; cx < nzeta; cx++) {
      tmat[(cx + nfact)*numParam + nfact + rx] = Ta[rx * nzeta + cx];
      tmat[(cx + nfact + ckoffset)*numParam + nfact + ckoffset + rx] = Tc[rx * nzeta + cx];
    }
  }

  double *dest = dmat + 2 * numParam * numParam;

  // It is probably possible to do this more efficiently than dgemm
  // since we know that we only care about the lower triangle.
  // I'm not sure whether this is worth optimizing. TODO

  char normal = 'n';
  char transpose = 't';
  double one = 1;
  double zero = 0;
  F77_CALL(dgemm)(&normal, &normal, &numParam, &numParam, &numParam,
		  &one, tmat, &numParam, dmat, &numParam, &zero, dest, &numParam);
  F77_CALL(dgemm)(&normal, &transpose, &numParam, &numParam, &numParam,
		  &one, dest, &numParam, tmat, &numParam, &zero, dmat, &numParam);

  {
    int row=0;
    int col=0;
    for (int dx=0; dx <= hsize; dx++) {
      out[numParam + dx] = dmat[numParam * col + row];
      if (row == col) {
	col=0; ++row;
      } else {
	++col;
      }
    }
  }

  Free(dmat);
}

static void
irt_rpf_mdim_nrm_dTheta(const double *spec, const double *param,
			const double *where, const double *dir,
			double *grad, double *hess)
{
  int numDims = spec[RPF_ISpecDims];
  int outcomes = spec[RPF_ISpecOutcomes];
  const double *aa = param;
  Eigen::VectorXd num(outcomes);
  Eigen::VectorXd ak(outcomes);
  double discr = dotprod(param, where, numDims);
  _nominal_rawprob2(spec, param, where, discr, ak.data(), num.data());

  double den = 0;
  for (int kx=0; kx < outcomes; kx++) {
    den += num[kx];
  }

  Eigen::VectorXd P(outcomes);
  for (int kx=0; kx < outcomes; kx++) {
    P[kx] = num[kx]/den;
  }

  for(int jx=0; jx < numDims; jx++) {
	  Eigen::VectorXd jak(outcomes);
	  Eigen::VectorXd jak2(outcomes);
    for (int ax=0; ax < outcomes; ax++) {
      jak[ax] = ak[ax] * aa[jx];
      jak2[ax] = jak[ax] * jak[ax];
    }
    double numjak = dotprod(num.data(), jak.data(), outcomes);
    double numjakden2 = numjak / den;
    numjakden2 *= numjakden2;
    double numjak2den = dotprod(num.data(), jak2.data(), outcomes) / den;

    for(int ix=0; ix < outcomes; ix++) {
      grad[ix] += dir[jx] * (ak[ix] * aa[jx] * P[ix] - P[ix] * numjak / den);
      hess[ix] += dir[jx] * (ak[ix]*ak[ix] * aa[jx]*aa[jx] * P[ix] -
			     2 * ak[ix] * aa[jx] * P[ix] * numjak / den +
			     2 * P[ix] * numjakden2 - P[ix] * numjak2den);
    }
  }
}

static void
irt_rpf_mdim_nrm_rescale(const double *spec, double *param, const int *paramMask,
			 const double *mean, const double *cov)
{
  int numDims = spec[RPF_ISpecDims];
  int nzeta = spec[RPF_ISpecOutcomes] - 1;
  double *alpha = param + numDims;
  double *gamma = param + numDims + nzeta;
  const double *Ta  = spec + RPF_ISpecCount;
  const double *Tc  = spec + RPF_ISpecCount + nzeta * nzeta;
  const double *iTc = spec + RPF_ISpecCount + 3 * nzeta * nzeta;

  double madj = dotprod(param, mean, numDims);

  for (int d1=0; d1 < numDims; d1++) {
    if (paramMask[d1] < 0) continue;
    param[d1] = dotprod(param+d1, cov + d1 * numDims + d1, numDims-d1);
  }

  Eigen::VectorXd ak(nzeta);
  ak.setZero();
  Eigen::VectorXd ck(nzeta);
  ck.setZero();

  for (int kx=0; kx < nzeta; kx++) {
    for (int tx=0; tx < nzeta; tx++) {
      int Tcell = tx * nzeta + kx;
      ak[kx] += Ta[Tcell] * alpha[tx];
      ck[kx] += Tc[Tcell] * gamma[tx];
    }
  }

  for (int kx=0; kx < nzeta; kx++) {
    ck[kx] += madj * ak[kx];
  }

  for (int kx=0; kx < nzeta; kx++) {
    int px = numDims + nzeta + kx;
    if (paramMask[px] < 0) continue;

    param[px] = 0;

    for (int tx=0; tx < nzeta; tx++) {
      int Tcell = tx * nzeta + kx;
      param[px] += iTc[Tcell] * ck[tx];
    }
  }
}


/********************************************************************************/
// Experimenting with LMP model - by CFF

static int
irt_rpf_1dim_lmp_numSpec(const double *spec)
{ return RPF_ISpecCount; }

static int
irt_rpf_1dim_lmp_numParam(const double *spec)
{
  int k = spec[RPF_ISpecCount];
  return(2+2*k);
}

// not sure what this does
static void
irt_rpf_1dim_lmp_paramInfo(const double *spec, const int param,
			   const char **type, double *upper, double *lower)
{
        //int k = spec[RPF_ISpecCount];
	*upper = nan("unset");
	*lower = nan("unset");

	*type = NULL;
	if (param == 0) {
		*type = "omega";
		//*lower = 1e-6;
	} else if (param == 1) {
		*type = "intercept";
	} else if (param %2 == 0 ){
		*type = "alpha";
	} else {
		*type = "tau";
	}
}

// Original R version did for multiple theta or x values (mp.val.noxi)
// Computes the value of m(x) given coefficients "b"
// and value(s) for th, omitting the intercept
static void
_mp_val(const double *th, const double *b, const int k,  double *out)
{
  int order = 2*k+1;
  out[0] = 0;
  for(int i=0; i<order; i++){
    out[0]+=b[i]*pow(*th,i+1);
  }
}

// Converts "a" coefficients for m'(x) to "b" coefficients for m(x)
// Omitting intercept xi
static void
_mp_getb (const double *a, const int k, double *b)
{
  int order = 2*k+1;
  for(int i=0; i<order; i++){
    b[i] = a[i]/(i+1);
  }
}

// Given item parameters, compute coefficeints "a" for m'(x)
// Does not include special case of k=0
// k - controls order of polynomial
// omega - sort of like slope parameter
// alpha - vector of length k that has one set of parameters
// tau - vector of length k that has another set of parameters
// dalpha - use derivative of T w.r.t. alpha in computations (0=no;1=first order; 2=second order)
// dtau - use derivative of T w.r.t. tau in computations (0=no;1=first order;2=second order)
// a_{k-1} - vector of lower order polynomial coefficients
static void
_mp_geta (const int k, const double *alpha, const double *tau,
	  const int dalpha, const int dtau, const Eigen::VectorXd *a, Eigen::VectorXd *newa)
{
  (*newa).setZero();
  int i, j;
  int indx = 0;
  int indx2 = 0;

  double beta = exp(*tau);
  double *t;

  // Part of Tk; the rest is 0's
  // optionally use derivatives wrt item parameters
  if(dalpha>0 && dtau>0){
    t = new double[3]{0, 0, 0};
  } else if (dalpha==1){
    t = new double[3]{0,-2,2.0*(*alpha)};
  } else if (dalpha==2){
    t = new double[3]{0,0,2};
  } else if (dtau==1||dtau==2){
    t = new double[3]{0,0,beta};
  } else {
    t = new double[3]{1, -2.0* (*alpha), pow(*alpha,2.0) + beta};
  }
  // a_k = T_k %*% a_{k-1} w/o saving T_k
  for(i=0; i<2*k-1; i++){
    for(j=0; j<2*k+1; j++){
       // only multiply by non-zero entries (unless derivatives are used, of course)
      if(j>=indx && j<indx+3){
        (*newa)[j]+=(*a)[i]*t[indx2];
	indx2++;
      }
    }
    indx++;
    indx2=0;
  }
}

// Recursively compute a starting at k=0 until desired k
// Optionally use derivatives in computation for da/deta where eta is an item parameter
// i.e., dalpha and dtau should be vectors with values 0-2 that determine
// which parameter derivatives are desired and what order of derivative
static void
_mp_getarec (const int k, const double *omega, const double *alpha, const double *tau,
	     int *dalpha, int *dtau, double *a)
{
  int i,j;
  //double olda[1];//Eigen::VectorXd olda(1);
  Eigen::VectorXd olda(1);
  olda[0] = exp(*omega);
  for(i=1;i<=k;i++){
    Eigen::VectorXd newa(i*2+1);
    _mp_geta(i,&alpha[i-1],&tau[i-1],dalpha[i-1],dtau[i-1],&olda,&newa);

    olda=newa; // need something like this after k=1
    // This is not optimal
    for(j=0;j<2*i+1;j++){
      a[j]=newa[j];
    }
  }
}

static void
_mp_dmda (const int k, const double *th, double *dmda)
{
  for(int i=0;i<(2*k+1);i++){
    dmda[i] = (1.0/(i+1.0))*pow(*th,i+1);
  }
}

static void
irt_rpf_1dim_lmp_prob(const double *spec,
		      const double *param, const double *th,
		      double *out)
{
  int k = spec[RPF_ISpecCount];
  double omega = param[0];
  double xi = param[1];
  double alpha[k];
  double tau[k];
  for(int i = 0; i<k; i++){
    alpha[i] = param[i*2+2];
    tau[i] = param[i*2+3];
  }

  double a[2*k+1];
  double b[2*k+1];
  double athb = 0;

  int* dalpha = new int[k];
  int* dtau = new int[k];
  std::fill(dalpha, dalpha+k,0);
  std::fill(dtau, dtau+k,0);

  _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
  _mp_getb(a,k,b);
  _mp_val(th,b,k,&athb);
  athb+=xi;

  if (athb < -EXP_STABLE_DOMAIN) athb = -EXP_STABLE_DOMAIN;
  else if (athb > EXP_STABLE_DOMAIN) athb = EXP_STABLE_DOMAIN;
  double pp = 1 / (1 + exp(-athb));
  out[0] = 1-pp;
  out[1] = pp;
}

static void
irt_rpf_1dim_lmp_deriv1(const double *spec,
				  const double *param,
				  const double *where,
				  const double *weight, double *out)
{
  int i, j;
  const int k = spec[RPF_ISpecCount];
  const int ord = 2*k+1;
  const int indxParam = 2+2*k;
  double omega = param[0];
  //double xi = param[1];
  double alpha[k];
  double tau[k];
  for(i = 0; i<k; i++){
    alpha[i] = param[i*2+2];
    tau[i] = param[i*2+3];
  }

  double QP[2];
  irt_rpf_1dim_lmp_prob(spec, param, where, QP);
  const double r1 = weight[1];
  const double r0 = weight[0];
  const double r1yP = QP[0]*r1;
  const double r0yP = -QP[1]*r0;
  const double PQ = QP[0]*QP[1];
  const double r1PQ = r1*PQ;
  const double r0PQ = r0*PQ;
  const double r1yP_r0yP = r1yP + r0yP;
  const double r1PQ_r0PQ = r1PQ + r0PQ;

  double dmdomega;
  double dmdalpha1;
  double dmdalpha2;
  double d2mdalpha1dalpha2;
  double dmdtau1;
  double dmdtau2;
  double d2mdtau1dtau2;
  double d2mdtaudalpha;
  double d2md2omega;
  double d2md2alpha;
  double d2md2tau;
  int* dalpha = new int[k];
  int* dtau = new int[k];
  std::fill(dalpha, dalpha+k,0);
  std::fill(dtau, dtau+k,0);

  double* a = new double[ord];
  //std::fill(a, a+ord, 0);

  double dmda[ord];
  _mp_dmda(k,where,dmda);

  // dldxi
  out[1] = r1yP_r0yP;

  // d2ld2xi
  out[hessianIndex(indxParam,1,1)] = -r1PQ_r0PQ;

  // dldomega
  _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
  dmdomega = dotprod(dmda, a, ord);
  d2md2omega = dmdomega; // same thing, at least with this parameterization
  out[0] = r1yP_r0yP*dmdomega;

  //d2ldxidomega
  out[hessianIndex(indxParam,1,0)] = -r1PQ_r0PQ*dmdomega;

  //d2ld2omega
  out[hessianIndex(indxParam,0,0)] = r1yP_r0yP*d2md2omega - r1PQ_r0PQ*dmdomega*dmdomega;


  for(i = 0; i<k; i++){
    // dldalpha
    dalpha[i] = 1;
    _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
    dalpha[i]=0;
    dmdalpha1 = dotprod(dmda, a, ord);
    out[i*2+2] = r1yP_r0yP*dmdalpha1;

    //d2ldalphadomega
    out[hessianIndex(indxParam,i*2+2,0)] = r1yP_r0yP*dmdalpha1-r1PQ_r0PQ*dmdalpha1*dmdomega;

    //d2ldxidalpha
    out[hessianIndex(indxParam,i*2+2,1)] = -r1PQ_r0PQ*dmdalpha1;

    // d2ld2alpha
    dalpha[i] = 2;
    _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
    dalpha[i] = 0;
    d2md2alpha = dotprod(dmda, a, ord);
    out[hessianIndex(indxParam,i*2+2,i*2+2)] = r1yP_r0yP*d2md2alpha - r1PQ_r0PQ*dmdalpha1*dmdalpha1;

    for(j = i+1; j<k; j++){
      //d2ldalpha1dalpha2
      dalpha[i]=0;
      dalpha[j]=1;
      _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
      dmdalpha2 = dotprod(dmda, a, ord);

      dalpha[i]=1;
      _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
      d2mdalpha1dalpha2 = dotprod(dmda, a, ord);

      out[hessianIndex(indxParam,j*2+2,i*2+2)] = r1yP_r0yP*d2mdalpha1dalpha2 - r1PQ_r0PQ*dmdalpha1*dmdalpha2;

      dalpha[i]=0;
      dalpha[j]=0;

    }

    // d2ldalpha1dalphatau
    for(j=0; j<k; j++){
      dalpha[i] = 0;
      dtau[j] = 1;
      _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
      dmdtau1 = dotprod(dmda, a, ord);

      dalpha[i]=1;
      _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
      d2mdtaudalpha = dotprod(dmda, a, ord);

      if(j>=i){
	out[hessianIndex(indxParam,j*2+3,i*2+2)] = r1yP_r0yP*d2mdtaudalpha-r1PQ_r0PQ*dmdalpha1*dmdtau1;
      } else {
	out[hessianIndex(indxParam,i*2+2,j*2+3)] = r1yP_r0yP*d2mdtaudalpha-r1PQ_r0PQ*dmdalpha1*dmdtau1;
      }

      dtau[j]=0;
      dalpha[i]=0;
    }
    //dalpha[i] = 0;

  }

  // dldtau
  for(i = 0; i<k; i++){
    dtau[i] = 1;
    _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
    dtau[i] = 0;
    dmdtau1 = dotprod(dmda, a, ord);
    out[i*2+3] = r1yP_r0yP*dmdtau1;

    //d2ldtaudomega
    out[hessianIndex(indxParam,i*2+3,0)] = r1yP_r0yP*dmdtau1-r1PQ_r0PQ*dmdtau1*dmdomega;

    //d2ldxidtau
    out[hessianIndex(indxParam,i*2+3,1)] = -r1PQ_r0PQ*dmdtau1;

    // d2ld2tau
    dtau[i] = 2;
    _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
    dtau[i] = 0;
    d2md2tau = dotprod(dmda, a, ord);
    out[hessianIndex(indxParam,i*2+3,i*2+3)] = r1yP_r0yP*d2md2tau - r1PQ_r0PQ*dmdtau1*dmdtau1;

    for(j = i+1; j<k; j++){
      //d2ldtau1dtau2
      dtau[i]=0;
      dtau[j]=1;
      _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
      dmdtau2 = dotprod(dmda, a, ord);

      dtau[i]=1;
      _mp_getarec(k, &omega, alpha, tau, dalpha, dtau, a);
      d2mdtau1dtau2 = dotprod(dmda, a, ord);

      out[hessianIndex(indxParam,j*2+3,i*2+3)] = r1yP_r0yP*d2mdtau1dtau2 - r1PQ_r0PQ*dmdtau1*dmdtau2;

      dtau[i]=0;
      dtau[j]=0;

    }
    //dtau[i] = 0;
  }
}

static void irt_rpf_1dim_lmp_deriv2(const double *spec,
				  const double *param,
				  double *out)
{

}


// Not edited yet
static void
irt_rpf_1dim_lmp_rescale(const double *spec, double *param, const int *paramMask,
			 const double *mean, const double *cov)
{
  /*double thresh = param[1] * -param[0];
  if (paramMask[0] >= 0) {
    param[0] *= cov[0];
  }
  if (paramMask[1] >= 0) {
    thresh += param[0] * mean[0];
    param[1] = thresh / -param[0];
    }*/
}
/********************************************************************************/

//static void noop() {}
static void notimplemented_deriv1(const double *spec,
				  const double *param,
				  const double *where,
				  const double *weight, double *out)
{ error("Not implemented"); }

static void notimplemented_deriv2(const double *spec,
				  const double *param,
				  double *out)
{ error("Not implemented"); }

static void notimplemented_dTheta(const double *spec, const double *param,
			const double *where, const double *dir,
			double *grad, double *hess)
{ error("Not implemented"); }


const struct rpf librpf_model[] = {
  { "drm1-",
    irt_rpf_1dim_drm_numSpec,
    irt_rpf_1dim_drm_numParam,
    irt_rpf_mdim_drm_paramInfo,
    irt_rpf_1dim_drm_prob,
    irt_rpf_logprob_adapter,
    notimplemented_deriv1,
    notimplemented_deriv2,
    irt_rpf_1dim_drm_dTheta,
    irt_rpf_1dim_drm_rescale,
  },
  { "drm",
    irt_rpf_mdim_drm_numSpec,
    irt_rpf_mdim_drm_numParam,
    irt_rpf_mdim_drm_paramInfo,
    irt_rpf_mdim_drm_prob,
    irt_rpf_logprob_adapter,
    irt_rpf_mdim_drm_deriv1,
    irt_rpf_mdim_drm_deriv2,
    irt_rpf_mdim_drm_dTheta,
    irt_rpf_mdim_drm_rescale,
  },
  { "grm",
    irt_rpf_mdim_grm_numSpec,
    irt_rpf_mdim_grm_numParam,
    irt_rpf_mdim_grm_paramInfo,
    irt_rpf_mdim_grm_prob,
    irt_rpf_logprob_adapter,
    irt_rpf_mdim_grm_deriv1,
    irt_rpf_mdim_grm_deriv2,
    irt_rpf_mdim_grm_dTheta,
    irt_rpf_mdim_grm_rescale,
  },
  { "nominal",
    irt_rpf_nominal_numSpec,
    irt_rpf_nominal_numParam,
    irt_rpf_nominal_paramInfo,
    irt_rpf_nominal_prob,
    irt_rpf_nominal_logprob,
    irt_rpf_nominal_deriv1,
    irt_rpf_nominal_deriv2,
    irt_rpf_mdim_nrm_dTheta,
    irt_rpf_mdim_nrm_rescale,
  },
  { "lmp",
    irt_rpf_1dim_lmp_numSpec,
    irt_rpf_1dim_lmp_numParam,
    irt_rpf_1dim_lmp_paramInfo,
    irt_rpf_1dim_lmp_prob,
    irt_rpf_logprob_adapter,
    irt_rpf_1dim_lmp_deriv1,
    irt_rpf_1dim_lmp_deriv2,  // not done yet
    notimplemented_dTheta, // not done yet
    irt_rpf_1dim_lmp_rescale, // not done yet
  }
};

const int librpf_numModels = (sizeof(librpf_model) / sizeof(struct rpf));

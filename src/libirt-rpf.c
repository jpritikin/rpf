/*
  Copyright 2012 Joshua Nathaniel Pritikin and contributors

  libirt-rpf is free software: you can redistribute it and/or modify
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

#include <math.h>
#include "libirt-rpf.h"

#ifndef M_LN2
#define M_LN2           0.693147180559945309417232121458        /* ln(2) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI    0.572364942924700087071713675677        /* log(sqrt(pi))
                                                                   == log(pi)/2 */
#endif

static const double DiscriminationSDLog = .5;
static const double LowerAsymptoteSD = .5;

/**
 * lognormal(x,sd) := (x*sd*sqrt(2*%pi))^-1 * exp(-(log(x)^2)/(2*sd^2))
 * log(lognormal(x,sd)), logexpand=super;
 */
static double lognormal_pdf(double aa)
{
	double sd2 = DiscriminationSDLog*DiscriminationSDLog;
	double loga = log(aa);
	return -loga*loga/(2*sd2) - loga - log(DiscriminationSDLog) - M_LN_SQRT_PI - M_LN2/2;
}

/**
 * diff(log(lognormal(a,sd))),a), logexpand=super;
 */
static double lognormal_gradient(double aa)
{
	double sd2 = DiscriminationSDLog*DiscriminationSDLog;
	return -log(aa)/(aa * sd2) - 1/aa;
}

static double logit(double prob)
{
	return log(prob/(1-prob));
}

// TODO Cai, Yang & Hansen (2011, p. 246)

/**
 * normal(x,mean,sd) := (sd*sqrt(2*%pi))^-1 * exp(-((x - mean)^2)/(2*sd^2))
 * log(normal(log(x/(1-x)),mean,sd)), logexpand=super;
 */
static double logitnormal_pdf(double cc, double mean, double sd)
{
	double sd2 = sd * sd;
	double offset = logit(cc) - mean;
	return -(offset*offset)/(2*sd2) - log(sd) - M_LN_SQRT_PI - M_LN2/2;
}

/**
 * factor(diff(log(normal(log(x/(1-x)),mean,sd)),x))
 */
static double logitnormal_gradient(double cc, double mean, double sd)
{
	double sd2 = sd * sd;
	return (logit(cc) - mean) / (sd2 * (cc-1) * cc);
}

int
irt_rpf_1dim_drm_numParam(const int numDims, const int numOutcomes)
{ return 3; }

void
irt_rpf_1dim_drm_logprob(const int numDims, const double *restrict param,
			 const double *restrict th,
			 const int numOutcomes, double *restrict out)
{
  double guessing = param[2];
  double pp = guessing + (1-guessing) / (1 + exp(-param[0] * (th[0] - param[1])));
  out[0] = log(1-pp);   // incorrect
  out[1] = log(pp);     // correct
}

double irt_rpf_1dim_drm_prior(const int numDims, const int numOutcomes,
			      const double *restrict param)
{
	double ll = lognormal_pdf(param[0]);
	double cc = param[2];
	if (cc > 0) {
		// configuration TODO
		ll += logitnormal_pdf(cc, logit(.25), LowerAsymptoteSD);
	}
	return ll;
}

/**
 * irf(a,b,c,th) := c+(1-c)/(1+exp(-a*(th - b)))
 * diff(log(1-irf(a,b,c,th)),a);  // 0
 * diff(log(irf(a,b,c,th)),a);    // 1
 * diff(log(1-irf(a,b,c,th)),b);  // 0
 * diff(log(irf(a,b,c,th)),b);    // 1
 * ratsimp(diff(log(1-irf(a,b,c,th)),c));
 * diff(log(irf(a,b,c,th)),c);
 */
void irt_rpf_1dim_drm_gradient(const int numDims, const int numOutcomes,
			       const double *restrict param, const int *paramMask,
			       const double *where, const double *weight, double *out)
{
	double aa = param[0];
	double bb = param[1];
	double cc = param[2];
	if (!where) {
		if (aa <= 0 || cc < 0 || cc >= 1) {
			out[0] = FP_NAN;
			out[1] = FP_NAN;
			out[2] = FP_NAN;
			return;
		}
		if (paramMask[0] >= 0) {
			out[0] *= (1-cc);
			out[0] += lognormal_gradient(aa);
		}
		if (paramMask[1] >= 0) {
			out[1] *= aa * (1-cc);
		}
		if (paramMask[2] >= 0) {
			if (cc == 0) { out[2] = FP_NAN; }
			else {
				out[2] += logitnormal_gradient(cc, logit(.25), LowerAsymptoteSD);
			}
		}
		return;
	}
	double th = where[0];
	double Eathb = exp(-aa * (th - bb));
	double Eathb2 = (Eathb+1)*(Eathb+1);
	double w0 = Eathb2 * (-(1-cc)/(Eathb+1) - cc + 1);
	double w1 = Eathb2 * ((1-cc)/(Eathb+1) + cc);
	if (paramMask[0] >= 0) {
		out[0] = (weight[0] * (bb-th)*Eathb / w0 +
			  weight[1] * -(bb-th)*Eathb / w1);
	}
	if (paramMask[1] >= 0) {
		out[1] = (weight[0] * Eathb / w0 +
			  weight[1] * -Eathb / w1);
	}
	if (paramMask[2] >= 0) {
		double ratio = (1-(1/(Eathb+1))) / ((1-cc)/(Eathb+1) + cc);
		out[2] = (weight[0] * (1/(cc-1)) +
			  weight[1] * ratio);
	}
}

int
irt_rpf_mdim_drm_numParam(const int numDims, const int numOutcomes)
{ return 2 + numDims; }

static double
dotprod(const double *v1, const double *v2, const int len)
{
  double dprod = 0;
  for (int dx=0; dx < len; dx++) {
    dprod += v1[dx] * v2[dx];
  }
  return dprod;
}

void
irt_rpf_mdim_drm_logprob(const int numDims, const double *restrict param,
			 const double *restrict th,
			 const int numOutcomes, double *restrict out)
{
  double dprod = dotprod(param, th, numDims);
  double diff = param[numDims];
  double guessing = param[numDims+1];
  double tmp = guessing + (1-guessing) / (1 + exp(-(dprod + diff)));
  out[0] = log(1-tmp);
  out[1] = log(tmp);
}

double irt_rpf_mdim_drm_prior(const int numDims, const int numOutcomes,
			      const double *restrict param)
{
	double ll=0;
	for (int dx=0; dx < numDims; dx++) {
		ll += lognormal_pdf(param[dx]);
	}
	double cc = param[numDims+1];
	if (cc > 0) {
		// configuration TODO
		ll += logitnormal_pdf(cc, logit(.25), LowerAsymptoteSD);
	}
	return ll;
}

/**
 * irf(a1,a2,b,c,th1,th2) := c+(1-c)/(1+exp(-(a1*th1+a2*th2+b)));
 * ratsimp(diff(log(1-irf(a1,a2,b,c,th1,th2)),a1)); //0
 * diff(log(1-irf(a1,a2,b,c,th1,th2)),a2);
 * diff(log(1-irf(a1,a2,b,c,th1,th2)),b);
 * diff(log(1-irf(a1,a2,b,c,th1,th2)),c);
 *
 * Wolfram Alpha
 * d/dx log(1- (c+(1-c)/(1+exp(-(x*f+y*g+z*h+b)))))
 * d/db log( (c+(1-c)/(1+exp(-(x*f+y*g+z*h+b)))))
 */
void irt_rpf_mdim_drm_gradient(const int numDims, const int numOutcomes,
			       const double *restrict param, const int *paramMask,
			       const double *where, const double *weight, double *out)
{
	const double *aa = param;
	double bb = param[numDims];
	double cc = param[numDims+1];
	if (!where) {
		if (cc < 0 || cc >= 1) {
			out[0] = FP_NAN;
			out[1] = FP_NAN;
			out[2] = FP_NAN;
			return;
		}
		for (int dx=0; dx < numDims; dx++) {
			if (aa[dx] <= 0) {
				out[0] = FP_NAN;
				out[1] = FP_NAN;
				out[2] = FP_NAN;
				return;
			}
		}
		for (int dx=0; dx < numDims; dx++) {
			if (paramMask[dx] >= 0) {
				out[dx] += lognormal_gradient(aa[dx]);
			}
		}
		if (paramMask[numDims] >= 0) {
			// OK
		}
		if (paramMask[numDims+1] >= 0) {
			if (cc == 0) { out[numDims+1] = FP_NAN; }
			else {
				out[numDims+1] += logitnormal_gradient(cc, logit(.25), LowerAsymptoteSD);
			}
		}
		return;
	}
	double Eathb = exp(dotprod(aa, where, numDims) + bb);
	for (int dx=0; dx < numDims; dx++) {
		if (paramMask[dx] >= 0) {
			out[dx] = (weight[0] * (where[dx]/(Eathb+1) - where[dx]) +
				   weight[1] * (where[dx]/(Eathb+1) - cc * where[dx] / (Eathb+cc)));
		}
	}
	if (paramMask[numDims] >= 0) {
		out[numDims] = (weight[0] * (1/(Eathb+1) - 1) +
				weight[1] * (-((cc-1)*Eathb)/((Eathb+1)*(Eathb+cc))));
	}
	if (paramMask[numDims+1] >= 0) {
		out[numDims+1] = (weight[0] * (1/(cc-1)) +
				  weight[1] * (1/(Eathb+cc)));
	}
}

int
irt_rpf_1dim_gpcm_numParam(const int numDims, const int numOutcomes)
{ return numOutcomes; }

void
irt_rpf_1dim_gpcm_logprob(const int numDims, const double *restrict param,
			  const double *restrict th,
			  const int numOutcomes, double *restrict out)
{
  double discr = param[0];
  double term1[numOutcomes];

  term1[numOutcomes - 1] = 0;

  for (int tx=0; tx < numOutcomes-1; tx++) {
    term1[tx] = -discr * (*th - param[tx+1]);
  }
  double sum = 0;
  double denom = 0;
  double term2[numOutcomes];
  for (int tx=numOutcomes-1; tx >= 0; tx--) {
    sum += term1[tx];
    term2[tx] = sum;
    denom += exp(sum);
  }
  denom = log(denom);
  for (int tx=0; tx < numOutcomes; tx++) {
    out[tx] = term2[tx] - denom;
  }
}

double irt_rpf_1dim_gpcm_prior(const int numDims, const int numOutcomes,
			       const double *restrict param)
{
	return lognormal_pdf(param[0]);
}

/**
 * Based on Muraki (1992, p. 167). I am not sure exactly how the
 * parameterization here relates to the parameterization in Muraki (1992).
 *
 * GPCM(3)
 * irf1(a,b1,b2,th) := 1/(exp(-a*(th-b2) + -a*(th-b1))+exp(-a*(th-b2))+ 1);
 * irf2(a,b1,b2,th) := (exp(-a*(th-b2)))/(exp(-a*(th-b2) + -a*(th-b1))+exp(-a*(th-b2))+ 1);
 * irf3(a,b1,b2,th) := (exp(-a*(th-b2) + -a*(th-b1)))/(exp(-a*(th-b2) + -a*(th-b1))+exp(-a*(th-b2))+ 1);
 * plot2d([irf1(1.5,-1,1,th), irf2(1.5,-1,1,th), irf3(1.5,-1,1,th)],[th,-3,3]);
 *
 * Muraki (1992) reparameterization?
 * P_0 irf1(a,b,d2,d1,th) := 1/(exp(-a*(th-b+d1) + -a*(th-b+d2))+exp(-a*(th-b+d1))+ 1)
 * P_1 irf2(a,b,d2,d1,th) := (exp(-a*(th-b+d1)))/(exp(-a*(th-b+d1) + -a*(th-b+d2))+exp(-a*(th-b+d1))+ 1)
 * P_2 irf3(a,b,d2,d1,th) := (exp(-a*(th-b+d1) + -a*(th-b+d2)))/(exp(-a*(th-b+d1) + -a*(th-b+d2))+exp(-a*(th-b+d1))+ 1)
 * plot2d([irf1(1.5,-1,0,-2,th), irf2(1.5,-1,0,-2,th), irf3(1.5,-1,0,-2,th)],[th,-3,3]);
 *
 * -b1 = -b+d2
 * -b2 = -b+d1
 */

static double
_1dim_gpcm_z(const int kk, const double *restrict param, double th)
{
  double sum=0;
  for (int vv=1; vv < kk; vv++) {
    sum += th - param[vv];
  }
  return param[0] * sum;
}

void irt_rpf_1dim_gpcm_gradient(const int numDims, const int numOutcomes,
				const double *restrict param, const int *paramMask,
				const double *where, const double *weight, double *out)
{
	double aa = param[0];
	if (!where) {
		if (aa <= 0) {
			out[0] = FP_NAN;
			out[1] = FP_NAN;
			out[2] = FP_NAN;
			return;
		}
		if (paramMask[0] >= 0) {
			out[0] /= aa;
			out[0] += lognormal_gradient(aa);
		}
		for (int bx=1; bx < numOutcomes; bx++) {
			if (paramMask[bx] >= 0) {
				out[bx] *= aa;
			}
		}
		return;
	}
	double pout[numOutcomes];
	irt_rpf_1dim_gpcm_logprob(numDims, param, where, numOutcomes, pout); // TODO add gpcm_prob
	for (int ox=0; ox < numOutcomes; ox++) { pout[ox] = exp(pout[ox]); }

	double th = where[0];

	if (paramMask[0] >= 0) {
		double grad=0;
		for (int kx=1; kx <= numOutcomes; kx++) {
			double sum=0;
			for (int cc=1; cc <= numOutcomes; cc++) {
				sum += _1dim_gpcm_z(cc, param, th) * pout[cc-1];
			}
			grad += weight[kx-1] * (_1dim_gpcm_z(kx, param, th) - sum);
		}
		out[0] = grad;
	}
	double total_weight=0;
	for (int cc=0; cc < numOutcomes; cc++) {
		total_weight += weight[cc];
	}
	for (int bx=1; bx < numOutcomes; bx++) {
		if (paramMask[bx] >= 0) {
			double grad = 0;
			for (int kx=0; kx <= bx-1; kx++) {
				grad += weight[kx] - pout[kx] * total_weight;
			}
			out[bx] = grad;
		}
	}
}

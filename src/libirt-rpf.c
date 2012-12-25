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

void
irt_rpf_1dim_drm_logprob(const double *restrict param,
			 const double th, double *restrict out)
{
  double guessing = param[2];
  out[1] = guessing + (1-guessing) / (1 + exp(-param[0] * (th - param[1])));
  out[0] = log(1-out[1]);
  out[1] = log(out[1]);
}

void
irt_rpf_mdim_drm_logprob(const int numDims, const double *restrict param,
			 const double *restrict th, double *restrict out)
{
  double dprod = 0;
  for (int dx=0; dx < numDims; dx++) {
    dprod += param[dx] * th[dx];
  }

  double diff = param[numDims];
  double guessing = param[numDims+1];
  double tmp = guessing + (1-guessing) / (1 + exp(-(dprod + diff)));
  out[0] = log(1-tmp);
  out[1] = log(tmp);
}

void
irt_rpf_1dim_gpcm_logprob(const int numOutcomes, const double *restrict param,
			  const double th, double *restrict out)
{
  double discr = param[0];
  double term1[numOutcomes];

  term1[numOutcomes - 1] = 0;

  for (int tx=0; tx < numOutcomes-1; tx++) {
    term1[tx] = -discr * (th - param[tx+1]);
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

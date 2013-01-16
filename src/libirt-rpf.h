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

#ifndef _LIBIRT_RPF_
#define _LIBIRT_RPF_

int
irt_rpf_1dim_drm_numParam(const int numDims, const int numOutcomes);
void
irt_rpf_1dim_drm_logprob(const int numDims, const double *restrict param,
			 const double *restrict th,
			 const int numOutcomes, double *restrict out);
double irt_rpf_1dim_drm_prior(const int numDims, const int numOutcomes,
			      const double *restrict param);
void irt_rpf_1dim_drm_gradient(const int numDims, const int numOutcomes,
			       const double *restrict param, const int *paramMask,
			       const double *where, const double *weight, double *out);

int
irt_rpf_mdim_drm_numParam(const int numDims, const int numOutcomes);
void
irt_rpf_mdim_drm_logprob(const int numDims, const double *restrict param,
			 const double *restrict th,
			 const int numOutcomes, double *restrict out);
double irt_rpf_mdim_drm_prior(const int numDims, const int numOutcomes,
			      const double *restrict param);
void irt_rpf_mdim_drm_gradient(const int numDims, const int numOutcomes,
			       const double *restrict param, const int *paramMask,
			       const double *where, const double *weight, double *out);

int
irt_rpf_1dim_gpcm_numParam(const int numDims, const int numOutcomes);
void
irt_rpf_1dim_gpcm_logprob(const int numDims, const double *restrict param,
			  const double *restrict th,
			  const int numOutcomes, double *restrict out);
double irt_rpf_1dim_gpcm_prior(const int numDims, const int numOutcomes,
			       const double *restrict param);
void irt_rpf_1dim_gpcm_gradient(const int numDims, const int numOutcomes,
				const double *restrict param, const int *paramMask,
				const double *where, const double *weight, double *out);

#endif

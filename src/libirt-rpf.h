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

void
irt_rpf_1dim_drm_logprob(const double *restrict param,
			 const double th, double *restrict out);

void
irt_rpf_mdim_drm_logprob(const int numDims, const double *restrict param,
			 const double *restrict th, double *restrict out);

void
irt_rpf_1dim_gpcm_logprob(const int numOutcomes, const double *restrict param,
			  const double th, double *restrict out);


#endif

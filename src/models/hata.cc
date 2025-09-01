/*****************************************************************************
*  HATA MODEL for Signal Server by Alex Farrant                              *
*  30 December 2013							     *
*  This program is free software; you can redistribute it and/or modify it   *
*  under the terms of the GNU General Public License as published by the     *
*  Free Software Foundation; either version 2 of the License or any later    *
*  version.								     *
*									     *
*  This program is distributed in the hope that it will useful, but WITHOUT  *
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     *
*  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License     *
*  for more details.							     *
*									     */

#include <cmath>

#include "hata.hh"

auto HATApathLoss(float f, float h_B, float h_M, float d, int mode) -> double
{
/*
HATA URBAN model for cellular planning
Frequency (MHz) 150 to 1500MHz
Base station height 30-200m
Mobile station height 1-10m
Distance 1-20km

mode 1 = URBAN
mode 2 = SUBURBAN
mode 3 = OPEN
*/
float lh_M = 0.0;
float C_H = 0.0;
const float logf = std::log10(f);

	if(f<200){
		lh_M = std::log10(1.54F * h_M);
		C_H = 8.29F * (lh_M * lh_M) - 1.1F;
	}else{
		lh_M = std::log10(11.75F * h_M);
		C_H = 3.2F * (lh_M * lh_M) - 4.97F;
	}

	const float L_u = 69.55F + 26.16F * logf - 13.82F * std::log10(h_B) - C_H + (44.9F - 6.55F * std::log10(h_B)) * std::log10(d);

	if (mode == 0 || mode == 1) {
		return L_u;	//URBAN
	}

	if (mode == 2) {	//SUBURBAN
		const float logf_28 = std::log10(f / 28);
		return L_u - 2 * logf_28 * logf_28 - 5.4;
	}

	if (mode == 3) {	//OPEN
		return L_u - 4.78 * logf * logf + 18.33 * logf - 40.94;
	}

	return 0;
}

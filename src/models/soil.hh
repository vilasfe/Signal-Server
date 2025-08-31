/*****************************************************************************
 *  Soil Path Loss model for Signal Server by Alex Farrant                    *
 *  21 February 2018                                                          *
 *                                                                            *
 *  This program is free software; you can redistribute it and/or modify it   *
 *  under the terms of the GNU General Public License as published by the     *
 *  Free Software Foundation; either version 2 of the License or any later    *
 *  version.                                                                  *
 *                                                                            *
 *  This program is distributed in the hope that it will useful, but WITHOUT  *
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     *
 *  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License     *
 *  for more details.                                                         *
 *
 * Frequency: Any MHz
 * Distance: Any Km
 * Terrain permittivity: 1 - 15 (Bad to Good)
 */

#ifndef _SOIL_HH_
#define _SOIL_HH_

#include "../common.h"

constexpr auto SoilPathLoss(float f, float d, float terdic) -> double
{
    const auto soil = (120/terdic);
    return(6.4 + _20log10(d) + _20log10(f)+(8.69*soil));
}

#endif /* _SOIL_HH_ */

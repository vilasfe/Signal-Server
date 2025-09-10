#ifndef _MAIN_HH_
#define _MAIN_HH_

#include <cmath>
#include <cstdio>
#include <string_view>

#include "common.h"

constexpr auto ReduceAngle(double angle) -> int
{
    /* This function normalizes the argument to
     *   an integer angle between 0 and 180 degrees */

    const double temp = std::acos(std::cos(angle * DEG2RAD));

    return static_cast<int>(std::rint(temp * RAD2DEG));
}

constexpr auto LonDiff(double lon1, double lon2) -> double
{
    /* This function returns the short path longitudinal
     *   difference between longitude1 and longitude2
     *   as an angle between -180.0 and +180.0 degrees.
     *   If lon1 is west of lon2, the result is positive.
     *   If lon1 is east of lon2, the result is negative. */

    double diff = lon1 - lon2;

    if (diff <= -180.0) {
        diff += 360.0;
    }

    if (diff >= 180.0) {
        diff -= 360.0;
    }

    return diff;
}

int PutMask(double lat, double lon, int value);
int OrMask(double lat, double lon, int value);
int GetMask(double lat, double lon);
void PutSignal(double lat, double lon, unsigned char signal);
unsigned char GetSignal(double lat, double lon);
double GetElevation(const struct site_t& location);
int AddElevation(double lat, double lon, double height, int size);

constexpr auto Distance(const struct site_t& site1, const struct site_t& site2) -> double
{
    /* This function returns the great circle distance
     *   in miles between any two site locations. */

    const double lat1 = site1.lat * DEG2RAD;
    const double lon1 = site1.lon * DEG2RAD;
    const double lat2 = site2.lat * DEG2RAD;
    const double lon2 = site2.lon * DEG2RAD;

    return
    3959.0 * std::acos(std::sin(lat1) * std::sin(lat2) +
    std::cos(lat1) * std::cos(lat2) * std::cos(lon1 - lon2));
}

constexpr auto Azimuth(const struct site_t& source, const struct site_t& destination) -> double {
    /* This function returns the azimuth (in degrees) to the
     destination as seen from* the location of the source. */

    const double dest_lat = destination.lat * DEG2RAD;
    const double dest_lon = destination.lon * DEG2RAD;

    const double src_lat = source.lat * DEG2RAD;
    const double src_lon = source.lon * DEG2RAD;

    const double sin_src_lat = std::sin(src_lat);
    const double cos_src_lat = std::cos(src_lat);
    const double sin_dest_lat = std::sin(dest_lat);

    /* Calculate Surface Distance */

    const double beta =
    std::acos(sin_src_lat * sin_dest_lat +
    cos_src_lat * std::cos(dest_lat) * std::cos(src_lon - dest_lon));

    /* Calculate Azimuth */

    const double num = sin_dest_lat - (sin_src_lat * std::cos(beta));
    const double den = cos_src_lat * std::sin(beta);
    /* Trap potential problems in acos() due to rounding */
    const double fraction = std::clamp(num / den, -1.0, 1.0);

    /* Calculate azimuth */

    double azimuth = std::acos(fraction);

    /* Reference it to True North */

    double diff = dest_lon - src_lon;

    if (diff <= -std::numbers::pi) {
        diff += TWOPI;
    }

    if (diff >= std::numbers::pi) {
        diff -= TWOPI;
    }

    if (diff > 0.0) {
        azimuth = TWOPI - azimuth;
    }

    return (azimuth * RAD2DEG);
}

double ElevationAngle(const struct site_t& source, const struct site_t& destination);
void ReadPath(const struct site_t& source, const struct site_t& destination);
double ElevationAngle2(const struct site_t& source, const struct site_t& destination, double er);
double ReadBearing(std::string_view input);
void ObstructionAnalysis(struct site_t xmtr, struct site_t rcvr, double f, FILE *outfile);

void free_elev();
void free_path();
void free_dem();
void alloc_elev();
void alloc_path();
void alloc_dem();
void do_allocs();

#endif /* _MAIN_HH_ */

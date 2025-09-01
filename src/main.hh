#ifndef _MAIN_HH_
#define _MAIN_HH_

#include <cmath>
#include <cstdio>

#include "common.h"

constexpr auto ReduceAngle(double angle) -> int
{
    /* This function normalizes the argument to
     *   an integer angle between 0 and 180 degrees */

    const double temp = std::acos(std::cos(angle * DEG2RAD));

    return static_cast<int>(std::rint(temp / DEG2RAD));
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

void *dec2dms(double decimal, char *string);
int PutMask(double lat, double lon, int value);
int OrMask(double lat, double lon, int value);
int GetMask(double lat, double lon);
void PutSignal(double lat, double lon, unsigned char signal);
unsigned char GetSignal(double lat, double lon);
double GetElevation(struct site location);
int AddElevation(double lat, double lon, double height, int size);
double Distance(struct site site1, struct site site2);
double Azimuth(struct site source, struct site destination);
double ElevationAngle(struct site source, struct site destination);
void ReadPath(struct site source, struct site destination);
double ElevationAngle2(struct site source, struct site destination, double er);
double ReadBearing(char *input);
void ObstructionAnalysis(struct site xmtr, struct site rcvr, double f, FILE *outfile);
void free_elev(void);
void free_path(void);
void free_dem(void);
void alloc_elev(void);
void alloc_path(void);
void alloc_dem(void);
void do_allocs(void);

#endif /* _MAIN_HH_ */

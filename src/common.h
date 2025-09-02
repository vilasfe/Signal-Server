#ifndef _COMMON_H_
#define _COMMON_H_

#include <cmath>
#include <numbers>
#include <string>

#define GAMMA 		2.5

constexpr double TWOPI = std::numbers::pi * 2;

constexpr double HALFPI = std::numbers::pi * 0.5;

constexpr double DEG2RAD = std::numbers::pi / 180.0;
constexpr double EARTHRADIUS = 20902230.97;
constexpr double METERS_PER_MILE = 1609.344;
constexpr double METERS_PER_FOOT = 0.3048;
constexpr double KM_PER_MILE = 1.609344;
constexpr double FEET_PER_MILE =5280.0;

constexpr double FOUR_THIRDS = 4.0/3.0;

struct dem {
	float min_north;
	float max_north;
	float min_west;
	float max_west;
	int max_el;
	int min_el;
	short **data;
	unsigned char **mask;
	unsigned char **signal;
};

struct site_t {
	double lat;
	double lon;
	float alt;
	char name[50];
	std::string filename;
};

struct path {
	double *lat;
	double *lon;
	double *elevation;
	double *distance;
	int length;
};

struct LR {
	double eps_dielect;
	double sgm_conductivity;
	double eno_ns_surfref;
	double frq_mhz;
	double conf;
	double rel;
	double erp;
	int radio_climate;
	int pol;
	float antenna_pattern[361][1001];
};

struct region {
	unsigned char color[128][3];
	int level[128];
	int levels;
};

extern int MAXPAGES;
extern int ARRAYSIZE;
extern int IPPD;

extern double min_north;
extern double max_north;
extern double min_west;
extern double max_west;
extern int ippd;
extern int MAXRAD;
extern int mpi;
extern int max_elevation;
extern int min_elevation;
extern int contour_threshold;
extern int loops;
extern int jgets;
extern int width;
extern int height;

extern double earthradius;
extern double north;
extern double east;
extern double south;
extern double west;
extern double max_range;
extern double dpp;
extern double ppd;
extern double yppd;
extern double fzone_clearance;
extern double clutter;
extern double dBm;
extern double loss;
extern double field_strength;
extern thread_local double *elev;
extern double westoffset;
extern double eastoffset;
extern double delta;
extern double cropLat;
extern double cropLon;

extern std::string sdf_path;
extern char gpsav;

extern unsigned char got_elevation_pattern;
extern unsigned char got_azimuth_pattern;
extern unsigned char metric;
extern unsigned char dbm;

extern struct dem *dem;
extern thread_local struct path path;
extern struct LR LR;
extern struct region region;

extern bool debug;

constexpr auto _10log10(auto x)
{
	return(4.342944f*std::log(x));
}

// use call with log/ln as this may be faster
// use constant of value 20.0/log(10.0)
constexpr auto _20log10(auto x)
{
	return(8.685889f*std::log(x));
}

#endif /* _COMMON_H_ */

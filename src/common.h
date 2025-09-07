#ifndef _COMMON_H_
#define _COMMON_H_

#include <cmath>
#include <numbers>
#include <string>

constexpr double GAMMA = 2.5;
constexpr double INV_GAMMA = 1 / GAMMA;

constexpr double TWOPI = std::numbers::pi * 2;

constexpr double HALFPI = std::numbers::pi * 0.5;

constexpr double DEG2RAD = std::numbers::pi / 180.0;
constexpr double RAD2DEG = 180.0 * std::numbers::inv_pi;
constexpr double EARTHRADIUS_FT = 20902230.97;
constexpr double METERS_PER_MILE = 1609.344;
constexpr double METERS_PER_FOOT = 0.3048;
constexpr double KM_PER_MILE = 1.609344;
constexpr double FEET_PER_MILE =5280.0;

constexpr double FOUR_THIRDS = 4.0/3.0;

struct dem_t {
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

struct path_t {
	double *lat;
	double *lon;
	double *elevation;
	double *distance;
	int length;
};

struct LR_t {
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

struct region_t {
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
extern int mpi;
extern int max_elevation;
extern int min_elevation;
extern int contour_threshold;
extern int width;
extern int height;

extern double earthradius;
extern double east;
extern double west;
extern double max_range;
extern double dpp;
extern double ppd;
extern double yppd;
extern double fzone_clearance;
extern double clutter;
extern thread_local double *elev;
extern double westoffset;
extern double eastoffset;
extern double delta;
extern double cropLat;
extern double cropLon;

extern std::string sdf_path;

extern unsigned char got_elevation_pattern;
extern unsigned char got_azimuth_pattern;
extern bool metric;
extern bool dbm;

extern struct dem_t *dem;
extern thread_local struct path_t path;
extern struct LR_t LR;
extern struct region_t region;

extern bool debug;

constexpr auto _10log10(auto&& x)
{
	return(4.342944F*std::log(x));
}

// use call with log/ln as this may be faster
// use constant of value 20.0/log(10.0)
constexpr auto _20log10(auto&& x)
{
	return(8.685889F*std::log(x));
}

#endif /* _COMMON_H_ */

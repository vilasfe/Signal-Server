#ifndef TILES_HH_
#define TILES_HH_

#include <string>
#include <string_view>

using tile_t = struct tile_t {
	std::string filename;
	int	cols = 0; // width
	int	rows = 0; // height
	double max_west = 0.0; // xll
	double min_north = 0.0; // yll
	double min_west = 0.0; // xur
	double max_north = 0.0; // yur
	double	cellsize = 0.0;
	long long datastart = 0LL;
	short	nodata = 0;
	short 	max_el = 0;
	short	min_el = 0;
	short	*data = nullptr;
	float 	precise_resolution = 0.0F;
	float	resolution = 0.0F;
	double	width_deg = 0.0;
	double	height_deg = 0.0;
	int		ppdx = 0;
	int		ppdy = 0;
};

auto tile_load_lidar(tile_t* tile, std::string_view filename) -> int;
auto tile_rescale(tile_t *tile, float scale) -> int;
void tile_destroy(tile_t *tile);

#endif

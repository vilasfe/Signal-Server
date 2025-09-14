#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <print>
#include <string>
#include <string_view>

#include "tiles.hh"

#include "common.h"

enum : std::uint16_t { MAX_LINE = 50000 };

auto tile_load_lidar(tile_t *tile, std::string_view filename) -> int {
	FILE *fd = nullptr;
	char line[MAX_LINE];
	char *pch = nullptr;

	/* Clear the tile data */
	*tile = tile_t{};

	/* Open the file handle and return on error */
	if ( fd = fopen(filename.data(),"r"); fd == nullptr ) {
		return errno;
	}

	/* This is where we read the header data */
	/* The string is split for readability but is parsed as a block */
	if( fscanf(fd,"%*s %d\n" "%*s %d\n" "%*s %lf\n" "%*s %lf\n" "%*s %lf\n" "%*s %d\n",&tile->cols,&tile->rows,&tile->max_west,&tile->min_north,&tile->cellsize,(int *)&tile->nodata) != 6 ){
		fclose(fd);
		return -1;
	}

	tile->datastart = ftell(fd);

	if(debug) {
		std::println(stderr, "w:{} h:{} s:{:f}", tile->cols, tile->rows, tile->cellsize);
	}

	/* Set the filename */
	tile->filename = filename;

	/* Perform xur calcs */
	tile->min_west = tile->max_west+(tile->cellsize*tile->cols);
	tile->max_north = tile->min_north+(tile->cellsize*tile->rows);

	eastoffset = std::max(eastoffset, tile->min_west);
	westoffset = std::min(westoffset, tile->max_west);

	 if (debug) {
	 	std::println(stderr,"{}, {} {:.7f}, {:.7f}, {:.7f}, {:.7f}, {:.7f}",tile->cols,tile->rows,tile->max_west,tile->min_north,tile->cellsize,tile->max_north,tile->min_west);
	 }

	// Greenwich straddling hack
	/* if (tile->xll <= 0 && tile->xur > 0) {
	 	tile->xll = (tile->xur - tile->xll); // full width
	 	tile->xur = 0.0; // budge it along so it's west of greenwich
	 	delta = eastoffset; // add to Tx longitude later
	 } else {*/
		// Transform WGS84 longitudes into 'west' values as society finishes east of Greenwich ;)
		if (tile->max_west >= 0) {
			tile->max_west = 360-tile->max_west;
		}
		if(tile->min_west >= 0) {
			tile->min_west = 360-tile->min_west;
		}
		if(tile->max_west < 0) {
			tile->max_west = -tile->max_west;
		}
		if(tile->min_west < 0) {
			tile->min_west = -tile->min_west;
		}
	// }

	if (debug) {
		std::println(stderr, "POST yll {:.7f} yur {:.7f} xur {:.7f} xll {:.7f} delta {:.6f}", tile->min_north, tile->max_north, tile->min_west, tile->max_west, delta);
	}

	/* Read the actual tile data */
	/* Allocate the array for the lidar data */
	if (tile->data = (short*) calloc(tile->cols * tile->rows, sizeof(short)); tile->data == nullptr ) {
		fclose(fd);
		tile->filename.clear();
		return ENOMEM;
	}

	size_t loaded = 0;
	for (size_t h = 0; h < static_cast<unsigned>(tile->rows); h++) {
		if (fgets(line, MAX_LINE, fd) != nullptr) {
			pch = strtok(line, " "); // split line into values
			for (size_t w = 0; w < static_cast<unsigned>(tile->cols) && pch != nullptr; w++) {
				/* If the data is less than a *magic* minimum, normalize it to zero */
				const auto nextval = static_cast<short>(std::max(0, std::atoi(pch)));
				tile->data[h*tile->cols + w] = nextval;
				loaded++;
				tile->max_el = std::max(tile->max_el, nextval);
				tile->min_el = std::min(tile->min_el, nextval);
				pch = strtok(nullptr, " ");
			}//while
		} else {
			std::println(stderr, "LIDAR error @ h {} file {}", h, filename);
		}//if
	}

	const double current_res_km = haversine_formula(tile->max_north, tile->max_west, tile->max_north, tile->min_west);
	tile->precise_resolution = static_cast<float>(current_res_km/std::max(tile->cols,tile->rows)*1000);

	// Round to nearest 0.5
	tile->resolution = tile->precise_resolution < 0.5F ? 0.5F : static_cast<float>(std::ceil((tile->precise_resolution * 2)+0.5) * 0.5);

	// Positive westing
	tile->width_deg = tile->max_west - tile->min_west >= 0 ? tile->max_west - tile->min_west : tile->max_west + (360 - tile->min_west);
	tile->height_deg = tile->max_north - tile->min_north;

	tile->ppdx = static_cast<int>(tile->cols / tile->width_deg);
	tile->ppdy = static_cast<int>(tile->rows / tile->height_deg);

	if (debug) {
		std::print(stderr,"Pixels loaded: {}/{} (PPD {}x{}, Res {:f} ({:.2f}))", loaded, tile->cols*tile->rows, tile->ppdx, tile->ppdy, tile->precise_resolution, tile->resolution);
	}

	/* All done, close the LIDAR file */
	fclose(fd);

	return 0;
}

/*
 * tile_rescale
 * This is used to resample tile data. It is particularly designed for
 * use with LIDAR tiles where the resolution can be anything up to 2m.
 * This function is capable of merging neighbouring pixel values
 * The scaling factor is the distance to merge pixels.
 * NOTE: This means that new resolutions can only increment in multiples of the original
 * (ie 2m LIDAR can be 4/6/8/... and 20m can be 40/60)
 */
auto tile_rescale(tile_t *tile, float scale) -> int {
	short *new_data = nullptr;
	size_t skip_count = 1;
	size_t copy_count = 1;

	if (scale == 1) {
		return 0;	
	}

	const size_t new_height = tile->rows * scale;
	const size_t new_width = tile->cols * scale;

	/* Allocate the array for the lidar data */
	if ( new_data = new short[new_height * new_width]; new_data == nullptr ) {
		return ENOMEM;
	}

	tile->max_el = std::numeric_limits<short>::min();
	tile->min_el = std::numeric_limits<short>::max();

	/* Making the tile data smaller */
	if (scale < 1) {
		skip_count = static_cast<size_t>(1 / scale);
	} else {
		copy_count = static_cast<size_t>(scale);
	}

	if (debug) {
		std::println(stderr,"Resampling tile {} [{:.1f}]:\n\tOld {}x{}. New {}x{}\n\tScale {:f} Skip {} Copy {}", tile->filename, tile->resolution, tile->cols, tile->rows, new_width, new_height, scale, skip_count, copy_count);
	}
	/* Nearest neighbour normalization. For each subsample of the original, simply
	 * assign the value in the top left to the new pixel 
	 * SOURCE: X / Y
	 * DEST:   I / J */

	for (size_t y = 0, j = 0; y < static_cast<unsigned>(tile->rows) && j < new_height; y += skip_count, j += copy_count) {

		for (size_t x = 0, i = 0; x < static_cast<unsigned>(tile->cols) && i < new_width; x += skip_count, i += copy_count) {
		
			/* These are for scaling up the data */
			for (size_t copy_y = 0; copy_y < copy_count; copy_y++) {
				for (size_t copy_x = 0; copy_x < copy_count; copy_x++) {
					const size_t new_j = j + copy_y;
					const size_t new_i = i + copy_x;
					/* Do the copy */
					new_data[ new_j * new_width + new_i ] = tile->data[y * tile->cols + x];
				}
			}
			/* Update local min / max values */
			tile->max_el = std::max(tile->max_el, tile->data[y * tile->cols + x]);
			tile->min_el = std::min(tile->min_el, tile->data[y * tile->cols + x]);
		}
	}

	/* Update the date in the tile */
	delete [] tile->data;
	tile->data = new_data;

	/* Update the height and width values */
	tile->rows = new_height;
	tile->cols = new_width;
	tile->resolution *= 1/scale;	// A scale of 2 is HALF the resolution
	tile->ppdx = tile->cols / tile->width_deg;
	tile->ppdy = tile->rows / tile->height_deg;
	// tile->width_deg *= scale;
	// tile->height_deg *= scale;
	if (debug) {
		std::println(stderr, "Resampling complete. New resolution: {:.1f}", tile->resolution);
	}

	return 0;
}

/*
 * tile_resize
 * This function works in conjuntion with resample_data. It takes a
 * resolution value in meters as its argument. It then calculates the
 * nearest (via averaging) resample value and calls resample_data
 */
auto tile_resize(tile_t* tile, int resolution) -> int {
	const double current_res_km = haversine_formula(tile->max_north, tile->max_west, tile->max_north, tile->min_west);
	const int current_res = static_cast<int>( std::ceil((current_res_km/IPPD)*1000) );
	const float scaling_factor = static_cast<float>(resolution) / current_res;
	if (debug) {
		std::println(stderr, "Resampling: Current {}m Desired {}m Scale {:.1f}", current_res, resolution, scaling_factor);
	}
	return tile_rescale(tile, scaling_factor);
}

/*
 * tile_destroy
 * This function simply destroys any data associated with a tile
 */
void tile_destroy(tile_t* tile) {
	delete [] tile->data;
}


#ifndef _TILES_HH_
#define _TILES_HH_

typedef struct _tile_t{
	char	*filename;
	union{
		int	cols = 0;
		int	width;
	};
	union{
		int	rows = 0;
		int	height;
	};
	union{
		double	xll = 0.0;
		double	max_west;
	};
	union{
		double	yll = 0.0;
		double	min_north;
	};
	union{
		double	xur = 0.0;
		double	min_west;
	};
	union{
		double	yur = 0.0;
		double	max_north;
	};
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
} tile_t, *ptile_t;

int tile_load_lidar(tile_t*, char *);
int tile_rescale(tile_t *, float);
void tile_destroy(tile_t *);

#endif

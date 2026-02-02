#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iosfwd>
#include <numbers>
#include <print>
#include <string>

#include "common.h"
#include "image.hh"
#include "inputs.hh"
#include "main.hh"
#include "outputs.hh"

#include "models/cost.hh"
#include "models/ecc33.hh"
#include "models/ericsson.hh"
#include "models/fspl.hh"
#include "models/hata.hh"
#include "models/itwom3.0.hh"
#include "models/sui.hh"

double Output::north = 0.0;
double Output::south = 0.0;
double Output::dBm = 0.0;
double Output::loss = 0.0;
double Output::field_strength = 0.0;
bool Output::gpsav = false;

int Output::height = 0;
int Output::width = 0;

void Output::DoPathLoss(std::string& filename, bool geo, bool kml, bool ngs, struct site_t *xmtr)
{
	/* This function generates a topographic map in Portable Pix Map
	   (PPM) format based on the content of flags held in the mask[][]
	   array (only).  The image created is rotated counter-clockwise
	   90 degrees from its representation in dem[][] so that north
	   points up and east points right in the image generated. */

	std::string mapfile;
	// Construction doesn't require passing the deleter again in C++20
	unique_file_ptr fd;

	auto ctx = Image::create(width, (kml ? height : height + 30), IMAGE_RGB, IMAGE_DEFAULT);
	int success = 0;

	const double conversion =
	    255.0 / std::pow(max_elevation - min_elevation,
			INV_GAMMA);

	if( success = Input::LoadLossColors(xmtr[0]); success != 0 ){
		std::println(stderr,"Error loading loss colors");
		exit(success);  // Now a fatal error!
	}

	if( !filename.empty() ) {

		if (filename[0] == 0) {
			filename = xmtr[0].filename;
			filename = filename.substr(0, filename.length() - 4);	/* Remove .qth */
		}

		mapfile = ctx->get_filename(std::string(filename));

		fd = unique_file_ptr(std::fopen(mapfile.data(), "wb"));

	} else {

		std::println(stderr,"Writing to stdout");
	}

	double minwest = min_west + dpp;

	if (minwest > 360.0) {
		minwest -= 360.0;
	}

	north = max_north - dpp;

	if (kml || geo) {
		south = min_north;	/* No bottom legend */
	} else {
		south = min_north - (30.0 / ppd);	/* 30 pixels for bottom legend */
	}

	east = (minwest < 180.0 ? -minwest : 360.0 - min_west);
	west = (max_west < 180 ? -max_west : 360 - max_west);

	if (debug) {
		std::println(stderr, "\nWriting \"{}\" ({}x{} pixmap image)...",
			!filename.empty() ? mapfile : "to stdout", width, (kml ? height : height + 30));
	}

	double lat = north;
	for (int y = 0; y < height; y++) {
		double lon = max_west;
		for (int x = 0; x < width; x++) {
			if (lon < 0.0) {
				lon += 360.0;
			}

			int x0 = 0;
			int y0 = 0;
			int indx = 0;
			bool found = false;
			for (; indx < MAXPAGES && !found;) {
				x0 = static_cast<int>(std::rint(ppd * (lat - static_cast<double>(dem[indx].min_north))));
				y0 = mpi - static_cast<int>(std::rint(ppd * (LonDiff(static_cast<double>(dem[indx].max_west), lon))));
				 // fix for multi-tile lidar
                              /*  if(width==10000 && (indx==1 || indx==3)){
                                        if(y0 >= 3432){ //3535
                                                y0=y0-3432;
                                        }
                                }*/

				if (x0 >= 0 && x0 <= mpi && y0 >= 0 && y0 <= mpi) {
					found = true;
				}
				else {
					indx++;
				}
			}

			if (found) {
				const unsigned char mask = dem[indx].mask[x0][y0];
				const int loss = dem[indx].signal[x0][y0];
				bool cityorcounty = false;

				int match = 255;

				unsigned red = 0;
				unsigned green = 0;
				unsigned blue = 0;

				if (loss <= region.level[0]) {
					match = 0;
				} else {
					for (int z = 1; (z < region.levels && match == 255); z++) {
						if (loss >= region.level[z - 1] && loss < region.level[z]) {
							match = z;
						}
					}
				}

				if (match < region.levels) {
					red = region.color[match][0];
					green = region.color[match][1];
					blue = region.color[match][2];
				}

				if ((mask & 2) != 0U) {
					/* Text Labels: Red or otherwise */

					if (red >= 180 && green <= 75 && blue <= 75 && loss == 0) {
						ctx->add_pixel(255 ^ red, 255 ^ green, 255 ^ blue);
					}
					else {
						ctx->add_pixel(255, 0, 0);
					}

					cityorcounty = true;
				}

				else if ((mask & 4) != 0U) {
					/* County Boundaries: Black */

					ctx->add_pixel(0, 0, 0);

					cityorcounty = true;
				}

				if (!cityorcounty) {
					if (loss == 0 || (contour_threshold != 0 && loss > std::abs(contour_threshold))) {
						if (ngs) {	/* No terrain */
							ctx->add_pixel(255, 255, 255);
						}
						else {
							/* Display land or sea elevation */

							if (dem[indx].data[x0][y0] == 0) {
								ctx->add_pixel(0, 0, 170);
							}
							else {
								auto terrain =
								    static_cast<unsigned>(
								    std::lround(
								     std::pow(static_cast<double>(dem[indx].data[x0][y0] - min_elevation), INV_GAMMA) * conversion));
								ctx->add_pixel(terrain, terrain, terrain);
							}
						}
					}

					else {
						/* Plot path loss in color */

						if (red != 0 || green != 0 || blue != 0) {
							ctx->add_pixel(red, green, blue);
						}
						else {	/* terrain / sea-level */

							if (dem[indx].data[x0][y0] == 0) {
								ctx->add_pixel(0, 0, 170);
							}
							else {
								/* Elevation: Greyscale */
								auto terrain =
								    static_cast<unsigned>(
								    std::lround(
								     std::pow(static_cast<double>(dem[indx].data[x0][y0] - min_elevation), INV_GAMMA) * conversion));
								ctx->add_pixel(terrain, terrain, terrain);
							}
						}
					}
				}
			}

			else {
				/* We should never get here, but if */
				/* we do, display the region as black */

				ctx->add_pixel(0, 0, 0);
			}
			lon = max_west - (dpp * (x+1));
		}
		lat = north - (dpp * (y+1));
	}

	if(success = ctx->write(!filename.empty() ? fd.get() : stdout); success != 0){
		std::println(stderr,"Error writing image");
		exit(success);
	}

}

auto Output::DoSigStr(std::string& filename, bool kml, bool ngs, struct site_t *xmtr) -> int
{
	/* This function generates a topographic map in Portable Pix Map
	   (PPM) format based on the signal strength values held in the
	   signal[][] array.  The image created is rotated counter-clockwise
	   90 degrees from its representation in dem[][] so that north
	   points up and east points right in the image generated. */

	std::string mapfile;
	FILE *fd = nullptr;
	auto ctx = Image::create(width, (kml ? height : height + 30), IMAGE_RGB, IMAGE_DEFAULT);
	int success = 0;

	const double conversion = 255.0 / std::pow(static_cast<double>(max_elevation - min_elevation), INV_GAMMA);

	if(success = Input::LoadSignalColors(xmtr[0]); success != 0 ){
		std::println(stderr,"Error loading signal colors");
		//exit(success);
	}

	if( !filename.empty() ) {
		if (filename[0] == 0) {
			filename = xmtr[0].filename;
			filename = filename.substr(0, filename.length() - 4);	/* Remove .qth */
		}

		mapfile = ctx->get_filename(filename);

		fd = fopen(mapfile.data(),"wb");
	} else {
		std::println(stderr,"Writing to stdout");
		fd = stdout;
	}

	double minwest = min_west + dpp;

	if (minwest > 360.0) {
		minwest -= 360.0;
	}

	north = max_north - dpp;

	south = min_north;	/* No bottom legend */

	east = (minwest < 180.0 ? -minwest : 360.0 - min_west);
	west = max_west < 180 ? -max_west : 360 - max_west;

	if (debug) {
		std::println(stderr, "\nWriting \"{}\" ({}x{} pixmap image)...",
			!filename.empty() ? mapfile : "to stdout", width, (kml ? height : height + 30));
	}

	double lat = north;
	for (int y = 0; y < height; y++) {
		double lon = max_west;
		for (int x = 0; x < width; x++) {
			if (lon < 0.0) {
				lon += 360.0;
			}

			int x0 = 0;
			int y0 = 0;
			bool found = false;
			int indx = 0;
			for (; indx < MAXPAGES && !found;) {
				x0 = static_cast<int>(std::rint(ppd *
					       (lat - static_cast<double>(dem[indx].min_north))));
				y0 = mpi -
				    static_cast<int>(std::rint(ppd *
					      (LonDiff(static_cast<double>(dem[indx].max_west), lon))));

				 // fix for multi-tile lidar
                           /*     if(width==10000 && (indx==1 || indx==3)){
                                        if(y0 >= 3432){ //3535
                                                y0=y0-3432;
                                        }
                                }
				*/

				if (x0 >= 0 && x0 <= mpi && y0 >= 0 && y0 <= mpi) {
					found = true;
				}
				else {
					indx++;
				}
			}

			if (found) {
				const unsigned char mask = dem[indx].mask[x0][y0];
				const int signal = (dem[indx].signal[x0][y0]) - 100;
				bool cityorcounty = false;
				int match = 255;

				unsigned red = 0;
				unsigned green = 0;
				unsigned blue = 0;

				if (signal >= region.level[0]) {
					match = 0;
				}
				else {
					for (int z = 1; (z < region.levels && match == 255); z++) {
						if (signal < region.level[z - 1] && signal >= region.level[z]) {
							match = z;
						}
					}
				}

				if (match < region.levels) {
					red = region.color[match][0];
					green = region.color[match][1];
					blue = region.color[match][2];
				}

				if ((mask & 2) != 0U) {
					/* Text Labels: Red or otherwise */

					if (red >= 180 && green <= 75 && blue <= 75) {
						ctx->add_pixel(255 ^ red, 255 ^ green, 255 ^ blue);
					}
					else {
						ctx->add_pixel(255, 0,	0);
					}

					cityorcounty = true;
				}

				else if ((mask & 4) != 0U) {
					/* County Boundaries: Black */

					ctx->add_pixel(0, 0, 0);

					cityorcounty = true;
				}

				if (!cityorcounty) {
					if (contour_threshold != 0 && signal < contour_threshold) {
						if (ngs) {
							ctx->add_pixel(255, 255, 255);
						}
						else {
							/* Display land or sea elevation */

							if (dem[indx].data[x0][y0] == 0) {
								ctx->add_pixel(0, 0, 170);
							}
							else {
								auto terrain =
								    static_cast<unsigned>(
								    std::lround(
								     std::pow(static_cast<double>(dem[indx].data[x0][y0] - min_elevation), INV_GAMMA) * conversion));
								ctx->add_pixel(terrain, terrain, terrain);
							}
						}
					}

					else {
						/* Plot field strength regions in color */

						if (red != 0 || green != 0 || blue != 0) {
							ctx->add_pixel(red, green, blue);
						}
						else {	/* terrain / sea-level */

							if (ngs) {
								ctx->add_pixel(255, 255, 255);
							}
							else {
								if (dem[indx].data[x0][y0] == 0) {
									ctx->add_pixel(0, 0, 170);
								}
								else {
									/* Elevation: Greyscale */
									auto terrain
									    =
									    static_cast<unsigned>(
									    std::lround(
									     std::pow
									     (static_cast<double>(dem[indx].data[x0][y0] - min_elevation), INV_GAMMA) * conversion));
									ctx->add_pixel(terrain, terrain, terrain);
								}
							}
						}
					}
				}
			}

			else {
				/* We should never get here, but if */
				/* we do, display the region as black */

				ctx->add_pixel(255, 255, 255);
			}
			lon = max_west - (dpp * (x+1));
		}
		lat = north - (dpp * (y+1));
	}

	if(success = ctx->write(fd); success != 0){
		std::println(stderr,"Error writing image");
		exit(success);
	}


	if( !filename.empty() ) {
		fclose(fd);
		fd = nullptr;
	}
	return 0;
}

void Output::DoRxdPwr(std::string filename, bool kml, bool ngs, struct site_t *xmtr)
{
	/* This function generates a topographic map in Portable Pix Map
	   (PPM) format based on the signal power level values held in the
	   signal[][] array.  The image created is rotated counter-clockwise
	   90 degrees from its representation in dem[][] so that north
	   points up and east points right in the image generated. */

	std::string mapfile;
	FILE *fd = stdout;
	auto ctx = Image::create(width, (kml ? height : height + 30), IMAGE_RGB, IMAGE_DEFAULT);
	int success = 0;

	const double conversion =
	    255.0 / std::pow(max_elevation - min_elevation,
			INV_GAMMA);

	if( success = Input::LoadDBMColors(xmtr[0]); success != 0 ){
		std::println(stderr,"Error loading DBM colors");
		exit(success);  //Now a fatal error!
	}

	if( !filename.empty() ) {

		if (filename[0] == 0) {
			filename = xmtr[0].filename;
			filename = filename.substr(0, filename.length() - 4);	/* Remove .qth */
		}

		mapfile = ctx->get_filename(filename);

		fd = fopen(mapfile.data(),"wb");

	} else {
		std::println(stderr,"Writing to stdout");
	}

	double minwest = min_west + dpp;

	if (minwest > 360.0) {
		minwest -= 360.0;
	}

	north = max_north - dpp;

	south = min_north;	/* No bottom legend */

	east = (minwest < 180.0 ? -minwest : 360.0 - min_west);
	west = (max_west < 180 ? -max_west : 360 - max_west);

	if (debug) {
		std::println(stderr, "\nWriting \"{}\" ({}x{} pixmap image)...",
			(!filename.empty() ? mapfile : "to stdout"), width, (kml ? height : height));
	}

	// Draw image of x by y pixels
	double lat = north;
	for (int y = 0; y < height; y++) {
		double lon = max_west;
		for (int x = 0; x < width; x++) {
			if (lon < 0.0) {
				lon += 360.0;
			}

			int x0 = 0;
			int y0 = 0;
			int indx = 0;
			bool found = false;
			for (; indx < MAXPAGES && !found;) {

				x0 = static_cast<int>(std::rint((ppd *
					      (lat - static_cast<double>(dem[indx].min_north)))));
				y0 = mpi -
				    static_cast<int>(std::rint(ppd *
					      (LonDiff(static_cast<double>(dem[indx].max_west),lon))));


				if (x0 >= 0 && x0 <= mpi && y0 >= 0 && y0 <= mpi) {
					found = true;
				}
				else {
					indx++;
				}


			}

			if (found) {
				const unsigned char mask = dem[indx].mask[x0][y0];
				const int dBm = (dem[indx].signal[x0][y0]) - 200;
				bool cityorcounty = false;
				int match = 255;

				unsigned red = 0;
				unsigned green = 0;
				unsigned blue = 0;

				if (dBm >= region.level[0]) {
					match = 0;
				}
				else {
					for (int z = 1; (z < region.levels && match == 255); z++) {
						if (dBm < region.level[z - 1] && dBm >= region.level[z]) {
							match = z;
						}
					}
				}

				if (match < region.levels) {
					red = region.color[match][0];
					green = region.color[match][1];
					blue = region.color[match][2];
				}

				if ((mask & 2) != 0U) {
					/* Text Labels: Red or otherwise */

					if (red >= 180 && green <= 75 && blue <= 75 && dBm != 0) {
						ctx->add_pixel(255 ^ red, 255 ^ green, 255 ^ blue);
					}
					else {
						ctx->add_pixel(255, 0, 0);
					}

					cityorcounty = true;
				}

				else if ((mask & 4) != 0U) {
					/* County Boundaries: Black */
					ctx->add_pixel(0, 0, 0);
					cityorcounty = true;
				}

				if (!cityorcounty) {
					if (contour_threshold != 0 && dBm < contour_threshold) {
						if (ngs) {	/* No terrain */
							ctx->add_pixel(255, 255, 255);
						}
						else {
							/* Display land or sea elevation */

							if (dem[indx].data[x0][y0] == 0) {
								ctx->add_pixel(0, 0, 170);
							}
							else {
								auto terrain =
								    static_cast<unsigned>(
								    std::lround(
								     std::pow(static_cast<double>(dem[indx].data[x0][y0] - min_elevation), INV_GAMMA) * conversion));
								ctx->add_pixel(terrain, terrain, terrain);
							}
						}
					}

					else {
						/* Plot signal power level regions in color */

						if (red != 0 || green != 0 || blue != 0) {
							ctx->add_pixel(red, green, blue);
						}
						else {	/* terrain / sea-level */

							if (ngs) {
								ctx->add_pixel(255, 255, 255); // WHITE
							}
							else {
								if (dem[indx].data[x0][y0] == 0) {
									ctx->add_pixel(0, 0, 170); // BLUE
								}
								else {
									/* Elevation: Greyscale */
									auto terrain
									    =
									    static_cast<unsigned>(
									    std::lround(
									     std::pow
									     (static_cast<double>(dem[indx].data[x0][y0] - min_elevation), INV_GAMMA) * conversion));
									ctx->add_pixel(terrain, terrain, terrain);
								}
							}
						}
					}
				}
			}

			else {
				/* We should never get here, but if */
				/* we do, display the region as black */

				ctx->add_pixel(255, 255, 255);
			}
			lon = max_west - (dpp * (x+1));
		}
		lat = north - (dpp * (y+1));
	}

	if(success = ctx->write(fd); success != 0){
		std::println(stderr,"Error writing image");
		exit(success);
	}

	fflush(fd);

	if( !filename.empty() ) {
		fclose(fd);
		fd = nullptr;
	}

}

void Output::DoLOS(std::string& filename, bool kml, bool ngs, struct site_t *xmtr)
{
	/* This function generates a topographic map in Portable Pix Map
	   (PPM) format based on the signal power level values held in the
	   signal[][] array.  The image created is rotated counter-clockwise
	   90 degrees from its representation in dem[][] so that north
	   points up and east points right in the image generated. */

	std::string mapfile;
	FILE *fd = nullptr;
	auto ctx = Image::create(width, (kml ? height : height + 30), IMAGE_RGB, IMAGE_DEFAULT);
	int success = 0;

	const double conversion =
	    255.0 / std::pow(max_elevation - min_elevation,
			INV_GAMMA);

	if( !filename.empty() ){

		if (filename[0] == 0) {
			filename = xmtr[0].filename;
			filename = filename.substr(0, filename.length() - 4);	/* Remove .qth */
		}

		mapfile = ctx->get_filename(filename);

		fd = fopen(mapfile.data(),"wb");

	} else {
		
		std::println(stderr,"Writing to stdout");
		fd = stdout;

	}

	double minwest = min_west + dpp;

	if (minwest > 360.0) {
		minwest -= 360.0;
	}

	north = max_north - dpp;

	south = min_north;	/* No bottom legend */

	east = (minwest < 180.0 ? -minwest : 360.0 - min_west);
	west = (max_west < 180 ? -max_west : 360 - max_west);

	if (debug) {
		std::println(stderr, "\nWriting \"{}\" ({}x{} pixmap image)...",
			!filename.empty() ? mapfile : "to stdout", width, (kml ? height : height + 30));
	}

	double lat = north;
	for (int y = 0; y < height; y++) {
		double lon = max_west;
		for (int x = 0; x < width; x++) {
			if (lon < 0.0) {
				lon += 360.0;
			}

			int x0 = 0;
			int y0 = 0;
			int indx = 0;
			bool found = false;
			for (; indx < MAXPAGES && !found;) {
				x0 = static_cast<int>(std::rint(ppd *
					       (lat - static_cast<double>(dem[indx].min_north))));
				y0 = mpi -
				    static_cast<int>(std::rint(ppd *
					      (LonDiff(static_cast<double>(dem[indx].max_west), lon))));

				if (x0 >= 0 && x0 <= mpi && y0 >= 0 && y0 <= mpi) {
					found = true;
				}
				else {
					indx++;
				}
			}

			if (found) {
				const unsigned char mask = dem[indx].mask[x0][y0];

				if ((mask & 2) != 0U) {
					/* Text Labels: Red */
					ctx->add_pixel(255, 0, 0);
				}
				else if ((mask & 4) != 0U) {
					/* County Boundaries: Light Cyan */
					ctx->add_pixel(128, 128, 255);
				}
				else {
					switch (mask & 57) {
					case 1:
						/* TX1: Green */
						ctx->add_pixel(0, 255, 0);
						break;

					case 8:
						/* TX2: Cyan */
						ctx->add_pixel(0, 255, 255);
						break;

					case 9:
						/* TX1 + TX2: Yellow */
						ctx->add_pixel(255, 255, 0);
						break;

					case 16:
						/* TX3: Medium Violet */
						ctx->add_pixel(147, 112, 219);
						break;

					case 17:
						/* TX1 + TX3: Pink */
						ctx->add_pixel(255, 192, 203);
						break;

					case 24:
						/* TX2 + TX3: Orange */
						ctx->add_pixel(255, 165, 0);
						break;

					case 25:
						/* TX1 + TX2 + TX3: Dark Green */
						ctx->add_pixel(0, 100, 0);
						break;

					case 32:
						/* TX4: Sienna 1 */
						ctx->add_pixel(255, 130, 71);
						break;

					case 33:
						/* TX1 + TX4: Green Yellow */
						ctx->add_pixel(173, 255, 47);
						break;

					case 40:
						/* TX2 + TX4: Dark Sea Green 1 */
						ctx->add_pixel(193, 255, 193);
						break;

					case 41:
						/* TX1 + TX2 + TX4: Blanched Almond */
						ctx->add_pixel(255, 235, 205);
						break;

					case 48:
						/* TX3 + TX4: Dark Turquoise */
						ctx->add_pixel(0, 206, 209);
						break;

					case 49:
						/* TX1 + TX3 + TX4: Medium Spring Green */
						ctx->add_pixel(0, 250, 154);
						break;

					case 56:
						/* TX2 + TX3 + TX4: Tan */
						ctx->add_pixel(210, 180, 140);
						break;

					case 57:
						/* TX1 + TX2 + TX3 + TX4: Gold2 */
						ctx->add_pixel(238, 201, 0);
						break;

					default:
						if (ngs) {	/* No terrain */
							ctx->add_pixel(255, 255, 255);
						}
						else {
							/* Sea-level: Medium Blue */
							if (dem[indx].data[x0][y0] == 0) {
								ctx->add_pixel(0, 0, 170);
							}
							else {
								/* Elevation: Greyscale */
								auto terrain =
								    static_cast<unsigned>(
								    std::lround(
								     std::pow(static_cast<double>(dem[indx].data[x0][y0] - min_elevation), INV_GAMMA) * conversion));
								ctx->add_pixel(terrain, terrain, terrain);
							}
						}
					}
				}
			}

			else {
				/* We should never get here, but if */
				/* we do, display the region as black */

				ctx->add_pixel(255, 255, 255);
			}
			lon = max_west - (dpp * (x+1));
		}
		lat = north - (dpp * (y+1));
	}

	if(success = ctx->write(fd); success != 0){
		std::println(stderr,"Error writing image");
		exit(success);
	}

	if( !filename.empty() ) {
		fclose(fd);
		fd = nullptr;
	}
}

auto Output::PathReport(struct site_t source, struct site_t destination, std::string& name,
		char graph_it, int propmodel, int pmenv, double rxGain) -> std::vector<double>
{
	/* This function writes a PPA Path Report (name.txt) to
	   the filesystem.  If (graph_it == 1), then gnuplot is invoked
	   to generate an appropriate output file indicating the Longley-Rice
	   model loss between the source and destination locations.
	   "filename" is the name assigned to the output file generated
	   by gnuplot.  The filename extension is used to set gnuplot's
	   terminal setting and output file type.  If no extension is
	   found, .png is assumed. */

	int errnum = 0;
	std::string basename;
	char term[30], ext[15],
	    report_name[80], block = 0;
	std::string strmode;
	double maxloss = -100000.0, minloss = 100000.0,
	    pattern = 1.0, patterndB = 0.0,
	    total_loss = 0.0, cos_xmtr_angle, cos_test_angle = 0.0,
	    source_alt, test_alt, dest_alt, source_alt2, dest_alt2,
	    distance, elevation,
	    free_space_loss = 0.0, eirp =
	    0.0, voltage, rxp, power_density, dkm;

	snprintf(report_name, 80, "%s.txt%c", name.data(), 0);
	const double four_thirds_earth = FOUR_THIRDS * EARTHRADIUS_FT;

	std::vector<double> elev(ARRAYSIZE + 10);

	auto fd2 = std::ofstream(report_name);

	std::println(fd2, "\n\t\t--==[ Path Profile Analysis ]==--\n");
	std::println(fd2, "Transmitter site: {}", source.name);

	if (source.lat >= 0.0) {

		if (source.lon <= 180){
			std::println(fd2, "Site location: {:.4f}, -{:.4f}",source.lat, source.lon);
		}else{
			std::println(fd2, "Site location: {:.4f}, {:.4f}",source.lat, 360 - source.lon);
		}
	}

	else {

		if (source.lon <= 180){
			std::println(fd2, "Site location: {:.4f}, -{:.4f}",source.lat, source.lon);
		}else{
			std::println(fd2, "Site location: {:.4f}, {:.4f}",source.lat, 360 - source.lon);
		}
	}

	if (metric) {
		std::println(fd2, "Ground elevation: {:.2f} meters AMSL",
			METERS_PER_FOOT * GetElevation(source));
		std::println(fd2,
			"Antenna height: {:.2f} meters AGL / {:.2f} meters AMSL",
			METERS_PER_FOOT * source.alt,
			METERS_PER_FOOT * (source.alt + GetElevation(source)));
	}

	else {
		std::println(fd2, "Ground elevation: {:.2f} feet AMSL", GetElevation(source));
		std::println(fd2, "Antenna height: {:.2f} feet AGL / {:.2f} feet AMSL",
			source.alt, source.alt + GetElevation(source));
	}

	double azimuth = Azimuth(source, destination);
	double angle1 = ElevationAngle(source, destination);
	double angle2 = ElevationAngle2(source, destination, earthradius);

	if (got_azimuth_pattern || got_elevation_pattern) {
		const int x = static_cast<int>(std::rint(10.0 * (10.0 - angle2)));

		if (x >= 0 && x <= 1000) {
			pattern =
			    static_cast<double>(LR.antenna_pattern[static_cast<int>(std::rint(azimuth))][x]);
		}

		patterndB = 20.0 * std::log10(pattern);
	}

	if (metric) {
		std::println(fd2, "Distance to {}: {:.2f} kilometers",
			destination.name, KM_PER_MILE * Distance(source, destination));
	}
	else {
		std::println(fd2, "Distance to {}: {:.2f} miles", destination.name, Distance(source, destination));
	}

	std::println(fd2, "Azimuth to {}: {:.2f} degrees grid", destination.name, azimuth);


	std::println(fd2, "Downtilt angle to {}: {:+.4f} degrees", destination.name, angle1);



	/* Receiver */

	std::println(fd2, "\nReceiver site: {}", destination.name);

	if (destination.lon <= 180){
		std::println(fd2, "Site location: {:.4f}, -{:.4f}",destination.lat, destination.lon);
	}else{
		std::println(fd2, "Site location: {:.4f}, {:.4f}",destination.lat, 360 - destination.lon);
	}

	if (metric) {
		std::println(fd2, "Ground elevation: {:.2f} meters AMSL",
			METERS_PER_FOOT * GetElevation(destination));
		std::println(fd2,
			"Antenna height: {:.2f} meters AGL / {:.2f} meters AMSL",
			METERS_PER_FOOT * destination.alt,
			METERS_PER_FOOT * (destination.alt +
					   GetElevation(destination)));
	}

	else {
		std::println(fd2, "Ground elevation: {:.2f} feet AMSL",
			GetElevation(destination));
		std::println(fd2, "Antenna height: {:.2f} feet AGL / {:.2f} feet AMSL",
			destination.alt,
			destination.alt + GetElevation(destination));
	}

	if (metric) {
		std::println(fd2, "Distance to {}: {:.2f} kilometers", source.name,
			KM_PER_MILE * Distance(source, destination));
	}
	else {
		std::println(fd2, "Distance to {}: {:.2f} miles", source.name,
			Distance(source, destination));
	}

	azimuth = Azimuth(destination, source);

	angle1 = ElevationAngle(destination, source);
	angle2 = ElevationAngle2(destination, source, earthradius);

	std::println(fd2, "Azimuth to {}: {:.2f} degrees grid", source.name, azimuth);


	std::println(fd2, "Downtilt angle to {}: {:+.4f} degrees", source.name, angle1);

	if (LR.frq_mhz > 0.0) {
		std::print(fd2, "\n\nPropagation model: ");

		switch (propmodel) {
		case 1:
			std::print(fd2, "Irregular Terrain Model");
			break;
		case 2:
			std::println(fd2, "Line of sight");
			break;
		case 3:
			std::println(fd2, "Okumura-Hata");
			break;
		case 4:
			std::println(fd2, "ECC33 (ITU-R P.529)");
			break;
		case 5:
			std::println(fd2, "Stanford University Interim");
			break;
		case 6:
			std::println(fd2, "COST231-Hata");
			break;
		case 7:
			std::println(fd2, "Free space path loss (ITU-R.525)");
			break;
		case 8:
			std::println(fd2, "ITWOM 3.0");
			break;
		case 9:
			std::println(fd2, "Ericsson");
			break;
		}

		std::print(fd2, "Model sub-type: ");

		switch (pmenv) {
		case 1:
			std::println(fd2, "City / Conservative");
			break;
		case 2:
			std::println(fd2, "Suburban / Average");
			break;
		case 3:
			std::println(fd2, "Rural / Optimistic");
			break;
		}
		std::println(fd2, "Earth's Dielectric Constant: {:.3f}", LR.eps_dielect);
		std::println(fd2, "Earth's Conductivity: {:.3f} Siemens/meter", LR.sgm_conductivity);
		std::println(fd2, "Atmospheric Bending Constant (N-units): {:.3f} ppm", LR.eno_ns_surfref);
		std::println(fd2, "Frequency: {:.3f} MHz", LR.frq_mhz);
		std::print(fd2, "Radio Climate: {} (", LR.radio_climate);

		switch (LR.radio_climate) {
		case 1:
			std::print(fd2, "Equatorial");
			break;

		case 2:
			std::print(fd2, "Continental Subtropical");
			break;

		case 3:
			std::print(fd2, "Maritime Subtropical");
			break;

		case 4:
			std::print(fd2, "Desert");
			break;

		case 5:
			std::print(fd2, "Continental Temperate");
			break;

		case 6:
			std::print(fd2, "Maritime Temperate, Over Land");
			break;

		case 7:
			std::print(fd2, "Maritime Temperate, Over Sea");
			break;

		default:
			std::print(fd2, "Unknown");
		}

		std::print(fd2, ")\nPolarisation: {} (", LR.pol);

		if (LR.pol == 0) {
			std::print(fd2, "Horizontal");
		}

		if (LR.pol == 1) {
			std::print(fd2, "Vertical");
		}

		std::println(fd2, ")\nFraction of Situations: {:.1f}{:c}", LR.conf * 100.0, 37);
		std::println(fd2, "Fraction of Time: {:.1f}{:c}", LR.rel * 100.0, 37);

		if (LR.erp != 0.0) {
			std::println(fd2, "\nReceiver gain: {:.1f} dBd / {:.1f} dBi", rxGain, rxGain+2.14);
			std::print(fd2, "Transmitter ERP plus Receiver gain: ");

			if (LR.erp < 1.0) {
				std::print(fd2, "{:.1f} milliwatts", 1000.0 * LR.erp);
			}

			if (LR.erp >= 1.0 && LR.erp < 10.0) {
				std::print(fd2, "{:.1f} Watts", LR.erp);
			}

			if (LR.erp >= 10.0 && LR.erp < 10.0e3) {
				std::print(fd2, "{:.0f} Watts", LR.erp);
			}

			if (LR.erp >= 10.0e3) {
				std::print(fd2, "{:.3f} kilowatts", LR.erp / 1.0e3);
			}

			dBm = 10.0 * (std::log10(LR.erp * 1000.0));
			std::println(fd2, " ({:+.2f} dBm)", dBm);
			std::println(fd2, "Transmitter ERP minus Receiver gain: {:.2f} dBm", dBm-rxGain);

			/* EIRP = ERP + 2.14 dB */

			std::print(fd2, "Transmitter EIRP plus Receiver gain: ");

			eirp = LR.erp * 1.636816521;

			if (eirp < 1.0) {
				std::print(fd2, "{:.1f} milliwatts", 1000.0 * eirp);
			}

			if (eirp >= 1.0 && eirp < 10.0) {
				std::print(fd2, "{:.1f} Watts", eirp);
			}

			if (eirp >= 10.0 && eirp < 10.0e3) {
				std::print(fd2, "{:.0f} Watts", eirp);
			}

			if (eirp >= 10.0e3) {
				std::print(fd2, "{:.3f} kilowatts", eirp / 1.0e3);
			}

			dBm = 10.0 * (std::log10(eirp * 1000.0));
			std::println(fd2, " ({:+.2f} dBm)", dBm);

			// Rx gain
			std::println(fd2, "Transmitter EIRP minus Receiver gain: {:.2f} dBm", dBm-rxGain);
		}

		std::println(fd2, "\nSummary for the link between {} and {}:\n", source.name, destination.name);

		if (patterndB != 0.0) {
			std::println(fd2, "{} antenna pattern towards {}: {:.3f} ({:.2f} dB)",
				source.name, destination.name, pattern, patterndB);
		}

		auto path = ReadPath(source, destination);	/* source=TX, destination=RX */

		/* Copy elevations plus clutter along
		   path into the elev[] array. */

		for (int x = 1; x < path.length - 1; x++) {
			elev[x + 2] =
			    METERS_PER_FOOT * (path.elevation[x] ==
					       0.0 ? path.
					       elevation[x] : (clutter +
							       path.
							       elevation[x]));
		}

		/* Copy ending points without clutter */

		elev[2] = path.elevation[0] * METERS_PER_FOOT;
		elev[path.length + 1] =
		    path.elevation[path.length - 1] * METERS_PER_FOOT;

		azimuth = std::rint(Azimuth(source, destination));

		for (int y = 2; y < (path.length - 1); y++) {	/* path.length-1 avoids LR error */
			distance = FEET_PER_MILE * path.distance[y];

			source_alt = four_thirds_earth + source.alt + path.elevation[0];
			dest_alt = four_thirds_earth + destination.alt +
			    path.elevation[y];
			dest_alt2 = dest_alt * dest_alt;
			source_alt2 = source_alt * source_alt;

			/* Calculate the cosine of the elevation of
			   the receiver as seen by the transmitter. */

			cos_xmtr_angle =
			    ((source_alt2) + (distance * distance) -
			     (dest_alt2)) / (2.0 * source_alt * distance);

			if (got_elevation_pattern) {
				/* If an antenna elevation pattern is available, the
				   following code determines the elevation angle to
				   the first obstruction along the path. */

				bool block = false;
				for (int x = 2; x < y && !block; x++) {
					distance =
					    FEET_PER_MILE * (path.distance[y] - path.distance[x]);
					test_alt =
					    four_thirds_earth +
					    path.elevation[x];

					/* Calculate the cosine of the elevation
					   angle of the terrain (test point)
					   as seen by the transmitter. */

					cos_test_angle =
					    ((source_alt2) +
					     (distance * distance) -
					     (test_alt * test_alt)) / (2.0 * source_alt * distance);

					/* Compare these two angles to determine if
					   an obstruction exists.  Since we're comparing
					   the cosines of these angles rather than
					   the angles themselves, the sense of the
					   following "if" statement is reversed from
					   what it would be if the angles themselves
					   were compared. */

					if (cos_xmtr_angle >= cos_test_angle) {
						block = true;
					}
				}

				/* At this point, we have the elevation angle
				   to the first obstruction (if it exists). */
			}

			//TODO: Seriously what is this array mess???
			/* Determine path loss for each point along the
			   path using Longley-Rice's point_to_point mode
			   starting at x=2 (number_of_points = 1), the
			   shortest distance terrain can play a role in
			   path loss. */

			elev[0] = y - 1;	/* (number of points - 1) */

			/* Distance between elevation samples */

			elev[1] =
			    METERS_PER_MILE * (path.distance[y] -
					       path.distance[y - 1]);

			/*
			   point_to_point(elev, source.alt*METERS_PER_FOOT,
			   destination.alt*METERS_PER_FOOT, LR.eps_dielect,
			   LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
			   LR.radio_climate, LR.pol, LR.conf, LR.rel, loss,
			   strmode, errnum);
			 */
			dkm = (elev[1] * elev[0]) / 1000;	// km

			switch (propmodel) {
			case 1:
				// Longley Rice ITM
				point_to_point_ITM(source.alt * METERS_PER_FOOT,
						   destination.alt *
						   METERS_PER_FOOT,
						   LR.eps_dielect,
						   LR.sgm_conductivity,
						   LR.eno_ns_surfref,
						   LR.frq_mhz, LR.radio_climate,
						   LR.pol, LR.conf, LR.rel, elev,
						   loss, strmode, errnum);
				break;
			case 3:
				//HATA 1, 2 & 3
				loss =
				    HATApathLoss(LR.frq_mhz, source.alt * METERS_PER_FOOT,
						 (path.elevation[y] * METERS_PER_FOOT) +
						 (destination.alt * METERS_PER_FOOT), dkm, pmenv);
				break;
			case 4:
				// COST231-HATA
				loss =
				    ECC33pathLoss(LR.frq_mhz, source.alt * METERS_PER_FOOT,
						  (path.elevation[y] * METERS_PER_FOOT) +
						  (destination.alt * METERS_PER_FOOT), dkm, pmenv);
				break;
			case 5:
				// SUI
				loss =
				    SUIpathLoss(LR.frq_mhz, source.alt * METERS_PER_FOOT,
						(path.elevation[y] * METERS_PER_FOOT) +
						(destination.alt * METERS_PER_FOOT), dkm, pmenv);
				break;
			case 6:
				loss =
				    COST231pathLoss(LR.frq_mhz, source.alt * METERS_PER_FOOT,
						    (path.elevation[y] * METERS_PER_FOOT) +
						    (destination.alt * METERS_PER_FOOT), dkm,pmenv);
				break;
			case 7:
				// ITU-R P.525 Free space path loss
				loss = FSPLpathLoss(LR.frq_mhz, dkm);
				break;
			case 8:
				// ITWOM 3.0
				point_to_point(source.alt * METERS_PER_FOOT,
					       destination.alt *
					       METERS_PER_FOOT, LR.eps_dielect,
					       LR.sgm_conductivity,
					       LR.eno_ns_surfref, LR.frq_mhz,
					       LR.radio_climate, LR.pol,
					       LR.conf, LR.rel, elev,
					       loss, strmode, errnum);
				break;
			case 9:
				// Ericsson
				loss =
				    EricssonpathLoss(LR.frq_mhz, source.alt * METERS_PER_FOOT,
						     (path.elevation[y] * METERS_PER_FOOT) +
						     (destination.alt *
						      METERS_PER_FOOT), dkm,
						     pmenv);
				break;

			default:
				point_to_point_ITM(source.alt * METERS_PER_FOOT,
						   destination.alt *
						   METERS_PER_FOOT,
						   LR.eps_dielect,
						   LR.sgm_conductivity,
						   LR.eno_ns_surfref,
						   LR.frq_mhz, LR.radio_climate,
						   LR.pol, LR.conf, LR.rel, elev,
						   loss, strmode, errnum);

			}

			if (!block) {
				elevation = (std::acos(cos_test_angle) / DEG2RAD) - 90.0;
			}
			else {
				elevation = (std::acos(cos_xmtr_angle) / DEG2RAD) - 90.0;
			}

			/* Integrate the antenna's radiation
			   pattern into the overall path loss. */

			const int x = static_cast<int>(std::rint(10.0 * (10.0 - elevation)));

			if (x >= 0 && x <= 1000) {
				pattern = static_cast<double>(LR.antenna_pattern[static_cast<int>(azimuth)][x]);

				if (pattern != 0.0){
					patterndB = 20.0 * std::log10(pattern);
				}else{
					patterndB = 0.0;
				}
			}

			else {
				patterndB = 0.0;
			}

			total_loss = loss - patterndB;

			maxloss = std::max(maxloss, total_loss);

			minloss = std::min(minloss, total_loss);

		}

		distance = Distance(source, destination);

		if (distance != 0.0) {
			free_space_loss =
			    36.6 + (20.0 * std::log10(LR.frq_mhz)) +
			    (20.0 * std::log10(distance));
			std::println(fd2, "Free space path loss: {:.2f} dB", free_space_loss);
		}

		std::println(fd2, "Computed path loss: {:.2f} dB", loss);


		if((loss*1.5) < free_space_loss) {
			std::println(fd2,"Model error! Computed loss of {:.1f}dB is greater than free space loss of {:.1f}dB. Check your inuts for model {}",loss,free_space_loss,propmodel);
			std::println(stderr,"Model error! Computed loss of {:.1f}dB is greater than free space loss of {:.1f}dB. Check your inuts for model {}",loss,free_space_loss,propmodel);
			return elev;
		}

		if (free_space_loss != 0.0) {
			std::println(fd2, "Attenuation due to terrain shielding: {:.2f} dB", loss - free_space_loss);
		}

		if (patterndB != 0.0) {
			std::println(fd2,"Total path loss including {} antenna pattern: {:.2f} dB", source.name, total_loss);
		}

		if (LR.erp != 0.0) {
			field_strength =
			    (139.4 + (20.0 * std::log10(LR.frq_mhz)) - total_loss) +
			    (10.0 * std::log10(LR.erp / 1000.0));

			/* dBm is referenced to EIRP */

			rxp = eirp / (std::pow(10.0, (total_loss / 10.0)));
			dBm = 10.0 * (std::log10(rxp * 1000.0));
			power_density =
			    (eirp /
			     (std::pow(10.0, (total_loss - free_space_loss) / 10.0)));
			/* divide by 4*PI*distance_in_meters squared */
			power_density /= (4.0 * std::numbers::pi * distance * distance * 2589988.11);

			std::println(fd2, "Field strength at {}: {:.2f} dBuV/meter", destination.name, field_strength);
			std::println(fd2, "Signal power level at {}: {:+.2f} dBm", destination.name, dBm);
			std::println(fd2, "Signal power density at {}: {:+.2f} dBW per square meter", destination.name, 10.0 * std::log10(power_density));
			voltage =
			    1.0e6 * std::sqrt(50.0 *
					 (eirp /
					  (std::pow(10.0, (total_loss - 2.14) / 10.0))));
			std::println(fd2, "Voltage across 50 ohm dipole at {}: {:.2f} uV ({:.2f} dBuV)",
				destination.name, voltage, 20.0 * std::log10(voltage));

			voltage =
			    1.0e6 * std::sqrt(75.0 *
					 (eirp /
					  (std::pow(10.0, (total_loss - 2.14) / 10.0))));
			std::println(fd2, "Voltage across 75 ohm dipole at {}: {:.2f} uV ({:.2f} dBuV)",
				destination.name, voltage, 20.0 * std::log10(voltage));
		}

		if (propmodel == 1) {
			std::print(fd2, "Longley-Rice model error number: {}",	errnum);

			switch (errnum) {
			case 0:
				std::println(fd2, " (No error)");
				break;

			case 1:
				std::println(fd2, "\n  Warning: Some parameters are nearly out of range.");
				std::println(fd2, "  Results should be used with caution.");
				break;

			case 2:
				std::println(fd2,
					"\n  Note: Default parameters have been substituted for impossible ones.");
				break;

			case 3:
				std::println(fd2,
					"\n  Warning: A combination of parameters is out of range for this model.");
				std::println(fd2,"  Results should be used with caution.");
				break;

			default:
				std::println(fd2,
					"\n  Warning: Some parameters are out of range for this model.");
				std::println(fd2,"  Results should be used with caution.");
			}
		}

	}

	ObstructionAnalysis(source, destination, LR.frq_mhz, fd2);
	fd2.close();

	std::print(stderr,
		"Path loss (dB), Received Power (dBm), Field strength (dBuV):\n{:.1f}\n{:.1f}\n{:.1f}",
		loss, dBm, field_strength);

	/* Skip plotting the graph if ONLY a path-loss report is needed. */

	if (graph_it) {
		if (name[0] == '.') {
			/* Default filename and output file type */

			basename = "profile";
			strncpy(term, "png\0", 4);
			strncpy(ext, "png\0", 4);
		}

		else {
			/* Extract extension and terminal type from "name" */

			ext[0] = 0;
			const int y = name.length();
			basename = name;

			int x = y - 1;
			for (; x > 0 && name[x] != '.'; x--) {}

			if (x > 0) {	/* Extension found */
				int z = x + 1;
				for (; z <= y && (z - (x + 1)) < 10; z++) {
					ext[z - (x + 1)] = tolower(name[z]);
					term[z - (x + 1)] = name[z];
				}

				ext[z - (x + 1)] = 0;	/* Ensure an ending 0 */
				term[z - (x + 1)] = 0;
				basename[x] = 0;
			}
		}

		if (ext[0] == 0) {	/* No extension -- Default is png */
			strncpy(term, "png\0", 4);
			strncpy(ext, "png\0", 4);
		}

		/* Either .ps or .postscript may be used
		   as an extension for postscript output. */

		if (strncmp(term, "postscript", 10) == 0) {
			strncpy(ext, "ps\0", 3);
		}
		else if (strncmp(ext, "ps", 2) == 0) {
			strncpy(term, "postscript enhanced color\0", 26);
		}

		auto fd = std::ofstream("ppa.gp");

		std::println(fd, "set grid");
		std::println(fd, "set yrange [{:2.3f} to {:2.3f}]", minloss, maxloss);
		std::println(fd, "set encoding iso_8859_1");
		std::println(fd, "set term {}", term);
		std::println(fd,
			"set title \"Path Loss Profile Along Path Between {} and {} ({:.2f}{:c} azimuth)\"",
			destination.name, source.name, Azimuth(destination,
							       source), 176);

		if (metric) {
			std::println(fd,
				"set xlabel \"Distance Between {} and {} ({:.2f} kilometers)\"",
				destination.name, source.name,
				KM_PER_MILE * Distance(destination, source));
		}
		else {
			std::println(fd,
				"set xlabel \"Distance Between {} and {} ({:.2f} miles)\"",
				destination.name, source.name,
				Distance(destination, source));
		}

		if (got_azimuth_pattern || got_elevation_pattern) {
			std::print(fd,
				"set ylabel \"Total Path Loss (including TX antenna pattern) (dB)");
		}
		else {
			std::print(fd, "set ylabel \"Longley-Rice Path Loss (dB)");
		}

		std::println(fd, "\"\nset output \"{}.{}\"", basename, ext);
		std::println(fd, "plot \"profile.gp\" title \"Path Loss\" with lines");

		fd.close();

		const int x = system("gnuplot ppa.gp");

		if (x != -1) {
			if (!gpsav) {
				//unlink("ppa.gp");
				//unlink("profile.gp");
				//unlink("reference.gp");
			}

		}

		else {
			std::println(stderr,
				"\n*** ERROR: Error occurred invoking gnuplot!");
		}
	}
	return elev;
}

void Output::SeriesData(const struct site_t& source, const struct site_t& destination, const std::string& name,
		bool fresnel_plot, bool normalised)
{
	double a, c, height = 0.0, cangle;
	double lambda = 0.0, f_zone =
	    0.0, fpt6_zone = 0.0, nm = 0.0, nb = 0.0, ed = 0.0, es = 0.0, r =
	    0.0, d = 0.0, d1 = 0.0, terrain;
	struct site_t remote;

	auto path = ReadPath(destination, source);
	const double azimuth = Azimuth(destination, source);
	const double distance = Distance(destination, source);
	const double refangle = ElevationAngle(destination, source);
	const double b = GetElevation(destination) + destination.alt + earthradius;

	if (debug) {
		std::println(stderr, "SeriesData: az = {:f}, dist = {:f}, ref = {:f}, b = {:f}", azimuth, distance, refangle, b);
	}
	
	if (fresnel_plot) {
		lambda = 9.8425e8 / (LR.frq_mhz * 1e6);
		d = FEET_PER_MILE * path.distance[path.length - 1];
	}

	if (normalised) {
		ed = GetElevation(destination);
		es = GetElevation(source);
		nb = -destination.alt - ed;
		nm = (-source.alt - es - nb) / (path.distance[path.length - 1]);
	}

	const std::string profilename = name + "_profile";
	const std::string referencename = name + "_reference";
	const std::string cluttername = name + "_clutter";
	const std::string curvaturename = name + "_curvature";
	const std::string fresnelname = name + "_fresnel";
	const std::string fresnel60name = name + "_fresnel60";

	auto fd = std::ofstream(profilename, std::ios::binary);
	std::ofstream fd1;
	if (clutter > 0.0) {
		fd.open(cluttername, std::ios::binary);
	}
	auto fd2 = std::ofstream(referencename, std::ios::binary);
	auto fd5 = std::ofstream(curvaturename, std::ios::binary);

	std::ofstream fd3;
	std::ofstream fd4;
	if ((LR.frq_mhz >= 20.0) && (LR.frq_mhz <= 100000.0) && fresnel_plot) {
		fd3.open(fresnelname, std::ios::binary);
		fd4.open(fresnel60name, std::ios::binary);
	}

	for (int x = 0; x < path.length - 1; x++) {
		remote.lat = path.lat[x];
		remote.lon = path.lon[x];
		remote.alt = 0.0;
		terrain = GetElevation(remote);
		if (x == 0) {
			terrain += destination.alt;	/* RX antenna spike */
		}

		a = terrain + earthradius;
		cangle = FEET_PER_MILE * Distance(destination, remote) / earthradius;
		c = b * std::sin(refangle * DEG2RAD + HALFPI) / std::sin(HALFPI -
							       refangle *
							       DEG2RAD -
							       cangle);
		height = a - c;

		/* Per Fink and Christiansen, Electronics
		 * Engineers' Handbook, 1989:
		 *
		 *   H = sqrt(lamba * d1 * (d - d1)/d)
		 *
		 * where H is the distance from the LOS
		 * path to the first Fresnel zone boundary.
		 */

		if ((LR.frq_mhz >= 20.0) && (LR.frq_mhz <= 100000.0) && fresnel_plot) {
			d1 = FEET_PER_MILE * path.distance[x];
			f_zone = -1.0 * std::sqrt(lambda * d1 * (d - d1) / d);
			fpt6_zone = f_zone * fzone_clearance;
		}

		if (normalised) {
			r = -(nm * path.distance[x]) - nb;
			height += r;

			if ((LR.frq_mhz >= 20.0) && (LR.frq_mhz <= 100000.0) && fresnel_plot) {
				f_zone += r;
				fpt6_zone += r;
			}
		}

		else {
			r = 0.0;
		}

		if (metric) {
			if (METERS_PER_FOOT * height > 0) {
				std::println(fd, "{:.3f} {:.3f}", KM_PER_MILE * path.distance[x], METERS_PER_FOOT * height);
			}

			if (fd1.is_open() && x > 0 && x < path.length - 2) {
				std::println(fd1, "{:.3f} {:.3f}", KM_PER_MILE * path.distance[x], METERS_PER_FOOT * (terrain == 0.0 ? height : (height + clutter)));
			}

			std::println(fd2, "{:.3f} {:.3f}", KM_PER_MILE * path.distance[x], METERS_PER_FOOT * r);
			std::println(fd5, "{:.3f} {:.3f}", KM_PER_MILE * path.distance[x], METERS_PER_FOOT * (height - terrain));
		}

		else {
			std::println(fd, "{:.3f} {:.3f}", path.distance[x], height);

			if (fd1.is_open() && x > 0 && x < path.length - 2) {
				std::println(fd1, "{:.3f} {:.3f}", path.distance[x],
					(terrain ==
					 0.0 ? height : (height + clutter)));
			}

			std::println(fd2, "{:.3f} {:.3f}", path.distance[x], r);
			std::println(fd5, "{:.3f} {:.3f}", path.distance[x], height - terrain);
		}

		if ((LR.frq_mhz >= 20.0) && (LR.frq_mhz <= 100000.0)
		    && fresnel_plot) {
			if (metric) {
				std::println(fd3, "{:.3f} {:.3f}",
					KM_PER_MILE * path.distance[x],
					METERS_PER_FOOT * f_zone);
				std::println(fd4, "{:.3f} {:.3f}",
					KM_PER_MILE * path.distance[x],
					METERS_PER_FOOT * fpt6_zone);
			}
			else {
				std::println(fd3, "{:.3f} {:.3f}", path.distance[x], f_zone);
				std::println(fd4, "{:.3f} {:.3f}", path.distance[x], fpt6_zone);
			}
		}
	}			// End of loop

	if (normalised) {
		r = -(nm * path.distance[path.length - 1]) - nb;
	}
	else {
		r = 0.0;
	}

	if (metric) {
		std::print(fd, "{:.3f} {:.3f}",
			KM_PER_MILE * path.distance[path.length - 1],
			METERS_PER_FOOT * r);
		std::print(fd2, "{:.3f} {:.3f}",
			KM_PER_MILE * path.distance[path.length - 1],
			METERS_PER_FOOT * r);
	}
	else {
		std::print(fd, "{:.3f} {:.3f}", path.distance[path.length - 1], r);
		std::print(fd2, "{:.3f} {:.3f}", path.distance[path.length - 1], r);
	}

	if ((LR.frq_mhz >= 20.0) && (LR.frq_mhz <= 100000.0) && fresnel_plot) {
		if (metric) {
			std::print(fd3, "{:.3f} {:.3f}",
				KM_PER_MILE * path.distance[path.length - 1],
				METERS_PER_FOOT * r);
			std::print(fd4, "{:.3f} {:.3f}",
				KM_PER_MILE * path.distance[path.length - 1],
				METERS_PER_FOOT * r);
		}
		else {
			std::print(fd3, "{:.3f} {:.3f}",
				path.distance[path.length - 1], r);
			std::print(fd4, "{:.3f} {:.3f}",
				path.distance[path.length - 1], r);
		}
	}

	std::println(stderr, "");
}

#ifndef _INPUTS_HH_
#define _INPUTS_HH_

#include <string_view>

#include <bzlib.h>
#include <zlib.h>

#include "common.h"

extern char scf_file[255];

class Input {
public:
    static int LoadSDF(char *name);
    static int LoadPAT(std::string_view az_filename, std::string_view el_filename);
    static int LoadSignalColors(struct site_t xmtr);
    static int LoadLossColors(struct site_t xmtr);
    static int LoadDBMColors(struct site_t xmtr);
    static int LoadTopoData(double max_lon, double min_lon, double max_lat, double min_lat);
    static int LoadUDT(std::string_view filename);
    static int loadLIDAR(const std::string& filenames, int resample);
    static int loadClutter(std::string_view filename, double radius, struct site_t tx);
    static int averageHeight(int h, int w, int x, int y);
    static constexpr char AZ_FILE_SUFFIX[] = ".az";
    static constexpr char EL_FILE_SUFFIX[] = ".el";
private:
    static int LoadSDF_SDF(char *name);
    static int LoadSDF_GZ(char *name);
    static int LoadSDF_BZ(char *name);
    static char *BZfgets(char *output, BZFILE *bzfd, unsigned length);
};

#endif /* _INPUTS_HH_ */

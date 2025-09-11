#ifndef INPUTS_HH_
#define INPUTS_HH_

#include <string>
#include <string_view>

#include <bzlib.h>
#include <zlib.h>

#include "common.h"

class Input {
public:
    static auto LoadSDF(char *name) -> int;
    static auto LoadPAT(std::string_view az_filename, std::string_view el_filename) -> int;
    static auto LoadSignalColors(struct site_t xmtr) -> int;
    static auto LoadLossColors(struct site_t xmtr) -> int;
    static auto LoadDBMColors(struct site_t xmtr) -> int;
    static auto LoadTopoData(double max_lon, double min_lon, double max_lat, double min_lat) -> int;
    static auto LoadUDT(std::string_view filename) -> int;
    static auto loadLIDAR(const std::string& filenames, int resample) -> int;
    static auto loadClutter(std::string_view filename, double radius, struct site_t tx) -> int;
    static auto averageHeight(int h, int w, int x, int y) -> int;
    static constexpr char AZ_FILE_SUFFIX[] = ".az";
    static constexpr char EL_FILE_SUFFIX[] = ".el";
    static int jgets;
    static std::string color_file;
private:
    static auto LoadSDF_SDF(char *name) -> int;
    static auto LoadSDF_GZ(char *name) -> int;
    static auto LoadSDF_BZ(char *name) -> int;
    static char *BZfgets(char *output, BZFILE *bzfd, unsigned length);
    static char *GZfgets(char *output, gzFile gzfd, unsigned length);
    static int bzerror;
    static int gzerr;
    static bool bzbuf_empty;
    static bool gzbuf_empty;
};

#endif /* INPUTS_HH_ */

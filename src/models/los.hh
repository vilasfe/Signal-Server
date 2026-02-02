#ifndef LOS_HH_
#define LOS_HH_

#include <memory>
#include <span>
#include <string>
#include <vector>

#include "../common.h"

class LOS {
public:
    static void PlotLOSMap(const struct site_t& source, double altitude, const std::string& plo_filename, bool use_threads);
    static void PlotPropagation(struct site_t source, double altitude, const std::string& plo_filename, int propmodel, int knifeedge, int haf, int pmenv, bool use_threads);
    static void PlotPath(const struct site_t& source, const struct site_t& destination, char mask_value);

    /*
     * Knife edge diffraction:
     * This is based upon a recognised formula like Huygens, but trades
     * thoroughness for increased speed which adds a proportional diffraction
     * effect to obstacles.
     */
    static constexpr auto ked(double freq, double rxh, double dkm, const std::span<double> elev) -> double
    {
        double rxobaoi = 0;
        double obh = 0;		// Obstacle height
        double obd = 0;		// Obstacle distance

        dkm = dkm * 1000;	// KM to metres

        // walk along path
        for (int n = 2; n < (dkm / elev[1]); n++) {

            const double d = (n - 2) * elev[1];	// no of points * delta = km

            //Find dip(s)
            if (elev[n] < obh) {
                // Angle from Rx point to obstacle
                rxobaoi = incidenceAngle((obh - (elev[n] + rxh)), d - obd);
            } else {
                // Line of sight or higher
                rxobaoi = 0;
            }

            //note the highest point
            if (elev[n] > obh) {
                obh = elev[n];
                obd = d;
            }
        }

        if (rxobaoi >= 0) {
            return (rxobaoi / (300 / freq))+3;	// Diffraction angle divided by wavelength (m)
        }
        return 1;
    }

private:
    static void PlotLOSPath(const struct site_t& source, const struct site_t& destination, unsigned char mask_value);
    [[nodiscard]] static auto PlotPropPath(const struct site_t& source, const struct site_t& destination, unsigned char mask_value, FILE * fd, int propmodel, int knifeedge, int pmenv) -> std::vector<double>;

    struct propagationRange {
        double min_west = 0.0;
        double max_west = 0.0;
        double min_north = 0.0;
        double max_north = 0.0;
        double altitude = 0.0;
        bool eastwest = false;
        bool los = true;
        bool use_threads = true;
        site_t source;
        unsigned char mask_value = 0;
        FILE *fd = nullptr;
        int propmodel = 0;
        int knifeedge = 0;
        int pmenv = 0;
    };

    static void rangePropagation(std::shared_ptr<propagationRange> v);

};

#endif /* LOS_HH_ */

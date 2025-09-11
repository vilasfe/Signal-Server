#ifndef LOS_HH_
#define LOS_HH_

#include <memory>
#include <string>

#include "../common.h"

class LOS {
public:
    static void PlotLOSMap(const struct site_t& source, double altitude, const std::string& plo_filename, bool use_threads);
    static void PlotPropagation(struct site_t source, double altitude, const std::string& plo_filename, int propmodel, int knifeedge, int haf, int pmenv, bool use_threads);
    static void PlotPath(const struct site_t& source, const struct site_t& destination, char mask_value);
private:
    static void PlotLOSPath(const struct site_t& source, const struct site_t& destination, unsigned char mask_value, FILE *fd);
    static void PlotPropPath(struct site_t source, struct site_t destination, unsigned char mask_value, FILE * fd, int propmodel, int knifeedge, int pmenv);

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

    static auto rangePropagation(std::shared_ptr<propagationRange> v) -> void*;
    static void beginThread(std::shared_ptr<propagationRange> arg);

    static auto ked(double freq, double rxh, double dkm) -> double;

    /*
     * Acute Angle from Rx point to an obstacle of height (opp) and
     * distance (adj)
     */
    static constexpr auto incidenceAngle(double opp, double adj) -> double
    {
        return std::atan2(opp, adj) * 180 * std::numbers::inv_pi;
    }

};

#endif /* LOS_HH_ */

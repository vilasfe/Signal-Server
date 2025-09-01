#ifndef _LOS_HH_
#define _LOS_HH_

#include <string>

#include "../common.h"

void PlotLOSPath(const struct site_t& source, const struct site_t& destination, unsigned char mask_value, FILE *fd);
void PlotPropPath(struct site_t source, struct site_t destination, unsigned char mask_value, FILE * fd, int propmodel, int knifeedge, int pmenv);
void PlotLOSMap(const struct site_t& source, double altitude, const std::string& plo_filename, bool use_threads);
void PlotPropagation(struct site_t source, double altitude, const std::string& plo_filename, int propmodel, int knifeedge, int haf, int pmenv, bool use_threads);
void PlotPath(const struct site_t& source, const struct site_t& destination, char mask_value);

#endif /* _LOS_HH_ */

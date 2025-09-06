#ifndef _OUTPUT_HH_
#define _OUTPUT_HH_

#include <string>

void DoPathLoss(std::string& filename, bool geo, bool kml,
		bool ngs, struct site_t *xmtr);
int DoSigStr(std::string& filename, bool geo, bool kml,
	      bool ngs, struct site_t *xmtr);
void DoRxdPwr(std::string filename, bool geo, bool kml,
	      bool ngs, struct site_t *xmtr);
void DoLOS(std::string& filename, bool geo, bool kml,
	   bool ngs, struct site_t *xmtr);
void PathReport(struct site_t source, struct site_t destination, std::string& name,
		char graph_it, int propmodel, int pmenv, double rxGain);
void SeriesData(struct site_t source, struct site_t destination, const std::string& name,
		bool fresnel_plot, bool normalised);

#endif /* _OUTPUT_HH_ */

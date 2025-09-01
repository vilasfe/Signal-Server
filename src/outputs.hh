#ifndef _OUTPUT_HH_
#define _OUTPUT_HH_

#include <string>

void DoPathLoss(std::string& filename, unsigned char geo, unsigned char kml,
		unsigned char ngs, struct site_t *xmtr, unsigned char txsites);
int DoSigStr(std::string& filename, unsigned char geo, unsigned char kml,
	      unsigned char ngs, struct site_t *xmtr, unsigned char txsites);
void DoRxdPwr(std::string filename, unsigned char geo, unsigned char kml,
	      unsigned char ngs, struct site_t *xmtr, unsigned char txsites);
void DoLOS(std::string& filename, unsigned char geo, unsigned char kml,
	   unsigned char ngs, struct site_t *xmtr, unsigned char txsites);
void PathReport(struct site_t source, struct site_t destination, std::string& name,
		char graph_it, int propmodel, int pmenv, double rxGain);
void SeriesData(struct site_t source, struct site_t destination, const std::string& name,
		unsigned char fresnel_plot, unsigned char normalised);

#endif /* _OUTPUT_HH_ */

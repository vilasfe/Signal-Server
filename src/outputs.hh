#ifndef _OUTPUT_HH_
#define _OUTPUT_HH_

#include <string>

class Output {
public:
	static void DoPathLoss(std::string& filename, bool geo, bool kml, bool ngs, struct site_t *xmtr);
	static int DoSigStr(std::string& filename, bool kml, bool ngs, struct site_t *xmtr);
	static void DoRxdPwr(std::string filename, bool kml, bool ngs, struct site_t *xmtr);
	static void DoLOS(std::string& filename, bool kml, bool ngs, struct site_t *xmtr);
	static void PathReport(struct site_t source, struct site_t destination, std::string& name,
		char graph_it, int propmodel, int pmenv, double rxGain);
	static void SeriesData(struct site_t source, struct site_t destination, const std::string& name,
		bool fresnel_plot, bool normalised);

private:
	static double north;
	static double south;
	static double dBm;
	static double loss;
	static double field_strength;
	static bool gpsav;
};

#endif /* _OUTPUT_HH_ */

#ifndef OUTPUT_HH_
#define OUTPUT_HH_

#include <span>
#include <string>
#include <vector>

class Output {
public:
	static void DoPathLoss(std::string& filename, bool geo, bool kml, bool ngs, struct site_t *xmtr);
	static auto DoSigStr(std::string& filename, bool kml, bool ngs, struct site_t *xmtr) -> int;
	static void DoRxdPwr(std::string filename, bool kml, bool ngs, struct site_t *xmtr);
	static void DoLOS(std::string& filename, bool kml, bool ngs, struct site_t *xmtr);
	[[nodiscard]] static auto PathReport(struct site_t source, struct site_t destination, std::string& name,
		char graph_it, int propmodel, int pmenv, double rxGain) -> std::vector<double>;
	static void SeriesData(const struct site_t& source, const struct site_t& destination, const std::string& name,
		bool fresnel_plot, bool normalised);

	static int width;
	static int height;

private:
	static double north;
	static double south;
	static double dBm;
	static double loss;
	static double field_strength;
	static bool gpsav;
};

#endif /* OUTPUT_HH_ */

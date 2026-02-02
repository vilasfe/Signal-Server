#ifndef ITWOM30_HH_
#define ITWOM30_HH_

#include <span>
#include <string>

#include "itm_types.hh"

class ITWOM3 {
public:

// TODO: return dbloss, make elev span<const double>, return string and errnum as std::expected or similar
static void point_to_point_ITM(double tht_m, double rht_m, double eps_dielect,
			double sgm_conductivity, double eno_ns_surfref,
			double frq_mhz, int radio_climate, int pol,
			double conf, double rel, std::span<double> elev,
			double &dbloss, std::string& strmode, int &errnum);
static void point_to_point(double tht_m, double rht_m, double eps_dielect,
		    double sgm_conductivity, double eno_ns_surfref,
		    double frq_mhz, int radio_climate, int pol, double conf,
		    double rel, std::span<double> elev,
			double &dbloss, std::string& strmode, int &errnum);

// Below this line public for unit testing access
[[nodiscard]] static auto lrprop(double d, prop_type & prop) -> propa_type;
[[nodiscard]] static auto lrprop2(double d, prop_type & prop) -> propa_type;

};

#endif /* ITWOM30_HH_ */

#ifndef ITWOM30_HH_
#define ITWOM30_HH_

#include <span>
#include <string>

void point_to_point_ITM(double tht_m, double rht_m, double eps_dielect,
			double sgm_conductivity, double eno_ns_surfref,
			double frq_mhz, int radio_climate, int pol,
			double conf, double rel, std::span<double> elev,
			double &dbloss, std::string& strmode, int &errnum);
void point_to_point(double tht_m, double rht_m, double eps_dielect,
		    double sgm_conductivity, double eno_ns_surfref,
		    double frq_mhz, int radio_climate, int pol, double conf,
		    double rel, std::span<double> elev,
			double &dbloss, std::string& strmode, int &errnum);

#endif /* ITWOM30_HH_ */

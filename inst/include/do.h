#ifndef DO_H
#define DO_H

#include <vector>

namespace NSM {

/**
	* Compute oxygen mass flux from concentration and discharge
	* @param DO dissolved oxygen concentration (mg/L)
	* @param Q discharge (m^3/min)
	* @return double; dissolved oxygen mass flux (mg/min)
*/
double oxyMassFlux(const std::vector<double> &DO, const std::vector<double> &Q);
double oxyMassFlux(double DO, double Q);

}

#endif
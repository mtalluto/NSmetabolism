#ifndef DO_H
#define DO_H

#include <vector>

namespace NSM {

/**
	* Compute oxygen mass flux from concentration and discharge
	* @param DO dissolved oxygen concentration (g/m^3)
	* @param Q discharge (m^3/s)
	* @return double; dissolved oxygen mass flux (g/min)
*/
double oxyMassFlux(const std::vector<double> &DO, const std::vector<double> &Q);
double oxyMassFlux(double DO, double Q);



/**
	* Compute advective concentration flux (g/[m^3 * min])
	* @param inputFlux Input DO mass flux (g/[m^3 * min])
 	* @param DOconc DO concentration (g/m^3)
	* @param Q discharge (m^3/s)
 	* @param area Cross sectional area (m)
 	* @param dx Segment length (m)
	* @return double; dissolved oxygen concentration flux (g/[m^3 * min])
*/
double advection(double inputFlux, double DOconc, double Q, double area, double dx);








#endif
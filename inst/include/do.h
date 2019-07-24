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


/**
	* Compute gross primary productivity from light (g O2 / [m^2 * day])
 	* @details P1 is in units of (W*day)/gO2, P2 is in (m^2 * day) / gO2
	* @param PAR Photosynthetically active radiation, in W/m^2 
	* @param lP1 The log of the slope of the PI curve
	* @param lP2 The log of the saturation term of the PI curve
 	* @references Uehlinger U. et al. 2000. Variability of photosynthesis-irradiance curves and
 	*  		ecosystem respiration in a small river. Freshwater Biology 44:493-507.
 	* @return double; computed GPP  
*/
double gpp(double PAR, double lP1, double lP2);
double gpp(double PAR, double lP1);


/**
	* Compute ecosystem respiration at in-situ temperature (gO2 / [m^2 * day])
	* @param temp Water temperature (degrees C)
	* @param ER24_20 Ecosystem respiration rate at 20 degrees (gO2 / [m^2 * day])
	* @return adjusted ER
*/
double er(double temp, double ER24_20);


/**
	* Computes reaeration flux (g/(m^2 * day))
	*
	* @param temp Water temperature (degrees C)
	* @param pressure Atmospheric pressure (hPa)
	* @param DO dissolved oxygen concentration, mg/L
	* @param k600 Gas transfer coefficient for Schmidt number of 600
	* @return computed rearation flux
	* @export
*/
double reaeration(double temp, double pressure, double DO, double k600);


/**
	* Compute gas transfer velocity for oxygen at a given temperature (m/day)
	* @param temp Water temperature (degrees C)
	* @param k600 Gas transfer coefficient for Schmidt number of 600 (m/day)
	* @references Wanninkhof R. (1992). Relationship between wind speed and gas exchange over the
	*    ocean. Journal of Geophysical Research, 97, 7373.\n
	*    Van de Bogert, M.C., Carpenter, S.R., Cole, J.J. & Pace, M.L. (2007). Assessing pelagic 
	*    and benthic metabolism using free water measurements. Limnology and Oceanography: 
	*    Methods, 5, 145–155.
	* @return k at the given temperature (m/day)
*/
double kT(double temp, double k600);

/**
	* Compute dissolved oxygen saturation given temperature and pressure (mg/L OR g/m^3)
	* @param temp Water temperature (degrees C)
	* @param P atmospheric pressure, in hPa
	* @references Benson BB and Krause D. 1984. The concentration and isotopic fractionation of
	*     oxygen dissolved in freshwater and seawater in equilibrium with the atmosphere. 
	*     Limnol. Oceanogr., 29, 620–632.
 	*     
 	* Benson BB and Krause D. 1980. The concentration and isotopic fractionation of gases
 	*     		dissolved in freshwater in equilibrium with the atmosphere. 1. Oxygen.  Limnol. 
 	*     		Oceanogr., 25(4):662-671
	* @return oxygen saturation concentration at given temperature and pressure
*/
double osat(double temp, double P);

}

#endif
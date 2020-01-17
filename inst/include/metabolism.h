#ifndef metabolism_h
#define metabolism_h

namespace NSM {

	/**
	 * @brief      Returns the log of the sum of exponentials
	 *
	 * @param	v1    First value to sum
	 * @param	v2    Second value to sum
	 *
	 * @return     { the log of the sum of the exponential of v1 and v2: log(exp(v1) + exp(v2)) }
	 */
	long double log_sum_exp(long double v1, long double v2);


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
	double gpp(double light, double lP1, double lP2);
	double gpp(double light, double lP1);

}

#endif
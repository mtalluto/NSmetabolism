/**
  * Computes reaeration flux
  *
  * @param temp Water temperature (degrees C)
  * @param pressure Atmospheric pressure (hPa)
  * @param DO dissolved oxygen concentration, mg/L
  * @param k600 Gas transfer coefficient for Schmidt number of 600
  * @return computed rearation flux
*/
real computeRF(real temp, real pressure, real DO, real k600) {
	return kT(temp, k600) * (osat(temp, pressure) - DO);
}


/**
    * Compute dissolved oxygen saturation given temperature and pressure (mg/L OR g/m^3)
    * @param temp Water temperature (degrees C)
    * @param P atmospheric pressure, in hPa
    * @references Benson BB and Krause D. 1984. The concentration and isotopic fractionation of
    *     oxygen dissolved in freshwater and seawater in equilibrium with the atmosphere. 
    *     Limnol. Oceanogr., 29, 620–632.
    *     
    * Benson BB and Krause D. 1980. The concentration and isotopic fractionation of gases
    *           dissolved in freshwater in equilibrium with the atmosphere. 1. Oxygen.  Limnol. 
    *           Oceanogr., 25(4):662-671
    * @return oxygen saturation concentration at given temperature and pressure
*/
real osat(real temp, real P) {  
	real Cstaro;
	real theta;
	real Pwv;
    real tempK = temp + 273.15;
	real hPaPerAtm = 1013.2500;
    real Patm = P / hPaPerAtm;

    // C*o, the unit standard atmospheric concentration by volume for oxygen; in mg/L
    // eqn 32 from Benson and Krause 1984.
	Cstaro = exp(-139.34411 + (1.575701e5 / tempK) - (6.642308e7 / tempK^2) + 
		  (1.243800e10 / tempK^3) - (8.621949e11 / tempK^4));

    // Benson and Krause 1980 eqn 13, 
    // the negative of the second pressure coefficient in the virial expansion for 
    // the real behavior of oxygen.
	theta = 0.000975 - 1.426e-5 * temp + 6.436e-8 * temp^2;

    // Benson and Krause 1980 eqn 23
    // the saturated vapor pressure of water in atmospheres at the temperature of equilibrium.
    // in atmospheres
	Pwv = exp(11.8571 - (3840.7 / tempK) - (216961 / tempK^2));

    // Benson and Krause 1980 eqn 28
	return Cstaro * ((Patm - Pwv) * (1 - theta * Patm)) / ((1 - Pwv) * (1 - theta));
}


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
real kT(real temp, real k600) {
	real Sc;
	// compute Schmidt number for oxygen
	// parameters from Wanninkhof 1992. appendix
	Sc = 1800.6 - 120.10 * temp + 3.7818 * temp^2 - 0.047608 * temp^3;
	
	// Van de Bogert et al eqn 5
	return k600 * (Sc / 600)^-0.5;
}


/**
  * Compute gross primary productivity from light (g O2 / [m^2 * day])
  * @details P1 is in units of (W*day)/gO2, P2 is in (m^2 * day) / gO2
  * @param lP1 The log of the slope of the PI curve
  * @param lP2 The log of the saturation term of the PI curve
  * @param data Fixed site data, including at least Q, area, and dx
  * @references Uehlinger U. et al. 2000. Variability of photosynthesis-irradiance curves and
  *         ecosystem respiration in a small river. Freshwater Biology 44:493-507.
  * @return double; computed GPP 
  * @param PAR Photosynthetically active radiation, in W/m^2 
  * @param lP1 The log of the slope of the PI curve
  * @param lP2 The log of the saturation term of the PI curve
*/
real computeGPP(real PAR, real lP1, real lP2) {
	real lGPP;
	if(PAR == 0)
		return 0;

	// Uehlinger et al 2000 eq 3b
	// GPP = PAR/(P1 + P2 * PAR)
	lGPP = log(PAR) - log_sum_exp(lP1, lP2 + log(PAR));

	return exp(lGPP);
}


/**
    * Compute ecosystem respiration at in-situ temperature (gO2 / [m^2 * day])
    * @param temp Water temperature (degrees C)
    * @param ER24_20 Ecosystem respiration rate at 20 degrees (gO2 / [m^2 * day])
    * @return adjusted ER
*/
real computeER(real temp, real ER24_20) {
	return ER24_20 * (1.045^(temp - 20));
}


/**
	* Compute expected value of k given slope and velocity, in m/day
	* @param slope Slope of the stream
	* @param velocity Velocity of the stream, m/s
	* @return real; expected value of k
*/
real k_mu(real slope, real velocity) {
	real kmu;
	
	kmu = 1162 * pow(slope, 0.77) * pow(velocity, 0.85); // in meters/day
	//kmu = kmu / (24*60); // convert to m/s 
	return kmu;
}

/**
	* Compute standard deviation of k given slope and velocity, in m/day
	* @param slope Slope of the stream
	* @param velocity Velocity of the stream, m/s
	* @return real; expected standard deviation of k
*/
real k_sd(real slope, real velocity) {
	real ksd;
	real a = 0.4595774;
	real b_s = 3.9006739;
	real b_v = 35.7944234;
	real b_sv = 270.8763474;
	
	ksd = a + b_s * slope + b_v * velocity + b_sv * slope * velocity;
	//ksd = ksd / (24*60); // convert to m/s 
	return ksd;
}

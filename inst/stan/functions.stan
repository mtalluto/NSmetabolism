functions{

/**
  * Compute oxygen saturation
  * @param temp Water temperature (degrees C)
  * @param P atmospheric pressure, in hPa
  * references Benson BB and Krause D. 1984. The concentration and isotopic fractionation of
  *     oxygen dissolved in freshwater and seawater in equilibrium with the atmosphere. 
  *     Limnol. Oceanogr., 29, 620–632.
  * @return oxygen saturation concentration at given temperature and pressure
*/
real osat(real temp, real P) {  
	real tempK;
	real Patm;
	real Cstaro;
	real theta;
	real Pwv;
	real hPaPerAtm = 1013.2500;

	tempK = temp + 273.15;
	Patm = P / hPaPerAtm;

	// C*o, the unit standard atmospheric concentration by volume for oxygen; in mg/kg
	// eqn 31 from Benson and Krause 1984.
	Cstaro = exp(-1.3874202e2 + (1.572288e5 / tempK) - (6.637149e7 / tempK^2) + 
		  (1.243678e10 / tempK^3) - (8.621061e11 / tempK^4));

	// eqn 13, 
	// the negative of the second pressure coefficient in the virial expansion for  gas 
	// the real behavior of oxygen.
	theta = 0.000975 - 1.426e-5 * temp + 6.436e-8 * temp^2;

	// eqn 23
	// the saturated vapor pressure of water in atmospheres at the temperature of equilibrium.
	// in atmospheres
	Pwv = exp(11.8571 - (3840.7 / tempK) - (216961 / tempK^2));

	// eqn 28
	return Cstaro * ((Patm - Pwv) * (1 - theta * Patm)) / ((1 - Pwv) * (1 - theta));
}


/**
  * Compute transport component
  * @param inputDO dissolved oxygen concentration of input water
  * @param outputDO dissolved oxygen concentration of focal pixel
  * @param Q discharge
  * @param area cross sectional area
  * @param dx length of pixel
*/
real computeAdvection(real inputDO, real outputDO, real Q, real area, real dx) {
	real inputMass = Q * inputDO;
	real outputMass = Q * outputDO;
	return (-1/area) * (outputMass - inputMass)/dx;
}



/**
  * Compute GPP from light
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
  * Compute ecosystem respiration at a given temperature
  * @param temp Water temperature (degrees C)
  * @param ER24_20 Ecosystem respiration rate for a 24-hour period at 20 degrees
*/
real computeER(real temp, real ER24_20) {
	return (ER24_20 / (60*24)) * (1.045^(temp - 20));
}

// }

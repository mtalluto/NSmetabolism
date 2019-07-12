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

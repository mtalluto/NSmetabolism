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
	* Compute expected value of k given slope and velocity, in m/s
	* @param slope Slope of the stream
	* @param velocity Velocity of the stream, m/s
	* @return real; expected value of k
*/
real k_mu(real slope, real velocity) {
	real kmu;
	
	kmu = 1162 * pow(slope, 0.77) * pow(velocity, 0.85); // in meters/day
	kmu = kmu / (24*60); // convert to m/s 
	return kmu;
}

/**
	* Compute standard deviation of k given slope and velocity, in m/s
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
	ksd = ksd / (24*60); // convert to m/s 
	return ksd;
}

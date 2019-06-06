real computeRF(real temp, real pressure, real DO, real k600) {
	return kT(temp, k600) * (osat(temp, pressure) - DO);
}

real osat(real temp, real P) {  
	real tempK;
	real Cstaro;
	real theta;
	real Pwv;

	tempK = temp + 273.15;

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
	return Cstaro * ((P - Pwv) * (1 - theta * P)) / ((1 - Pwv) * (1 - theta));
}

real computeAdvection(real inputDO, real outputDO, real Q, real area, real dx) {
	real inputMass = Q * inputDO;
	real outputMass = Q * outputDO;
	return (-1/area) * (outputMass - inputMass)/dx;
}

real kT(real temp, real k600) {
	real Sc;
	// compute Schmidt number for oxygen
	// parameters from Wanninkhof 1992. appendix
	Sc = 1800.6 - 120.10 * temp + 3.7818 * temp^2 - 0.047608 * temp^3;
	
	// Van de Bogert et al eqn 5
	return k600 * (Sc / 600)^-0.5;
}

real computeGPP(real PAR, real lP1, real lP2) {
	real lGPP;
	if(PAR == 0)
		return 0;

	// Uehlinger et al 2000 eq 3b
	// GPP = PAR/(P1 + P2 * PAR)
	lGPP = log(PAR) - log_sum_exp(lP1, lP2 + log(PAR));

	return exp(lGPP);
}

real computeER(real temp, real ER24_20) {
	return (ER24_20 / (60*24)) * (1.045^(temp - 20));
}
functions {
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
		if(lP2 == 0) {
			lGPP = log(PAR) - lP1;
		} else {
			lGPP = log(PAR) - log_sum_exp(lP1, lP2 + log(PAR));
		}
		return exp(lGPP);
	}

	real computeRF(real temp, real pressure, real DO, real k600) {
		return kT(temp, k600) * (osat(temp, pressure) - DO);
	}

	real computeER(real temp, real ER24_20) {
		return (ER24_20 / (60*24)) * (1.045^(temp - 20));
	}

	real computeAdvection(real inputDO, real outputDO, real discharge, real area, real dx) {
		real inputMass = discharge * inputDO;
		real outputMass = discharge * outputDO;
		return (-1/area) * (outputMass - inputMass)/dx;
	}

}
data {
	int<lower = 1> nDO; // DO number of observations;

	vector<lower=0> [nDO] DO; // DO Observations
	int<lower = 1> time [nDO]; // time of each DO observation, in minutes
	int<lower = 1, upper = 2> site [nDO]; // site of DO observations, 1 is upstream, 2 is down
	vector<lower=0> [2] DOinitial; // initial values for each site

	int<lower = max(time)> maxTime; // latest time at which we want predictions
	matrix [maxTime, 2] temp;
	matrix [maxTime, 2] PAR;

	// stream characteristics; should be average of two sites for scalars
	real <lower = 0> slope;
	real <lower = 0> velocity;
	real <lower = 0> pressure;
	real <lower = 0> discharge; // two station model assumes no lateral input, constant Q
	real <lower = 0> area; // cross sectional area
	real <lower = 0> dx;
	real <lower = 0> z;
	real <lower = 0> dt;
}
parameters {
	vector<lower=0> [2] lP1;
	vector<lower=0> [2] lP2; 
	vector<upper=0> [2] ER24_20;
	real<lower=0> k600;
	real<lower=0> sigma;
}
transformed parameters {
	matrix [maxTime,2] gpp = rep_matrix(0, maxTime, 2);
	matrix [maxTime,2] er = rep_matrix(0, maxTime, 2);
	matrix [maxTime,2] DO_pr; 

	DO_pr[1,1] = DOinitial[1];
	DO_pr[1,2] = DOinitial[2];
	for(i in 2:maxTime) {
		for(j in 1:2) {
			real ddodt;
			real rf;
			real inputDO;
			real adv;

			if(j == 1)
			// NOTE - input concentration for upstream site is just the concentration
			// of the upstream site from the time period before
				inputDO = DO_pr[i-1,j];
			else
				inputDO = DO_pr[i-1,j-1];
			adv = computeAdvection(inputDO, DO_pr[i-1,j], discharge, area, dx);
			rf = computeRF(temp[i-1,j], pressure, DO_pr[i-1,j], k600);
			gpp[i,j] = computeGPP(PAR[i-1,j], lP1[j], lP2[j]);
			er[i,j] = computeER(temp[i-1,j], ER24_20[j]);
			ddodt = adv + (gpp[i,j] + er[i,j] + rf) / z;
			DO_pr[i,j] = DO_pr[i-1,j] + ddodt * dt;
		}
	}
}
model {
	for(i in 1:nDO) {
		DO[i] ~ normal(DO_pr[time[i], site[i]], sigma);
	}

	k600 ~ normal((1162 * pow(slope, 0.77) * pow(velocity, 0.85))/(24*60), 0.0001462944 + 0.0012564499 * slope + 0.0124307051 * velocity + 0.0961094198 * slope * velocity);
	lP1 ~ normal(9, 1);
	ER24_20 ~ normal(0, 10);
}


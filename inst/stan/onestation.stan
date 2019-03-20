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

	real computeGPP(real PAR, real lP1) {
		// Uehlinger et al 2000 eq 3b
		return exp(log(PAR) - lP1);
	}
	// real computeGPP(real PAR, real lP1, real lP2) {
	// 	real lGPP;
	// 	if(PAR == 0)
	// 		return 0;

	// 	// Uehlinger et al 2000 eq 3b
	// 	// GPP = PAR/(P1 + P2 * PAR)
	// 	if(lP2 == 0) {
	// 		lGPP = log(PAR) - lP1;
	// 	} else {
	// 		lGPP = log(PAR) - log_sum_exp(lP1, lP2 + log(PAR));
	// 	}
	// 	return exp(lGPP);
	// }

	real computeRF(real temp, real pressure, real DO, real k600) {
		return kT(temp, k600) * (osat(temp, pressure) - DO);
	}

	real computeER(real temp, real ER24_20) {
		return (ER24_20 / (60*24)) * (1.045^(temp - 20));
	}
}
data {
	int<lower = 1> nDO; // DO number of observations;
	int<lower = 1> nTime; // number of integration time steps
	int<lower = 0> timesDO [nDO]; // time (minutes) of DO
	int<lower = 0> timesInt [nTime]; // integration times
	vector<lower = 0> [nTime] PAR;
	vector [nTime] temp;
	real<lower=0> z;
	real pressure;
	vector<lower=0> [nDO] DO;
	real slope;
	real velocity;
}
parameters {
	real<lower=0> lP1;
	// real<lower=0> lP2; 
	real<upper=0> ER24_20;
	real<lower=0> k600;
	real<lower=0> sigma;
}
transformed parameters {
	vector [nDO] gpp = rep_vector(0, nDO);
	vector [nDO] er = rep_vector(0, nDO);
	vector [nDO] doPredicted; 

	{
		real state;
		int DOindex = 2;
		state = DO[1];
		doPredicted[1] = DO[1];

		for(i in 2:nTime) {
			real rf = computeRF(temp[i-1], pressure, state, k600);
			real gppi = computeGPP(PAR[i-1], lP1);
			real eri = computeER(temp[i-1], ER24_20);
			real ddodt = (gppi + eri + rf)/z;
			real dt = timesInt[i] - timesInt[i - 1];

			state += ddodt * dt;
			gpp[DOindex] += gppi;
			er[DOindex] += eri;

			if(timesInt[i] == timesDO[DOindex]) {
				doPredicted[DOindex] = state;
				DOindex += 1;
			}
		}
	}
}
model {
	DO ~ normal(doPredicted, sigma);
	k600 ~ normal((1162 * pow(slope, 0.77) * pow(velocity, 0.85))/(24*60), 0.0001462944 + 0.0012564499 * slope + 0.0124307051 * velocity + 0.0961094198 * slope * velocity);
	lP1 ~ normal(9, 1);
	ER24_20 ~ normal(0, 10);
}


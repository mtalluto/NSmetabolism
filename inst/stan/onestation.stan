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
}
data {
	int<lower = 1> nDO; // DO number of observations;

	vector<lower=0> [nDO] DO; // DO Observations
	int<lower = 1> time [nDO]; // time of each DO observation, in minutes
	real<lower=0> DOinitial; // initial value DO

	int<lower = max(time)> maxTime; // latest time at which we want predictions
	vector [maxTime] temp;
	vector [maxTime] PAR;

	// stream characteristics; should be average of two sites for scalars
	real <lower = 0> slope;
	real <lower = 0> velocity;
	real <lower = 0> pressure;
	real <lower = 0> z; // depth
	real <lower = 0> dt; // length (in minutes) of a single time step
}

parameters {
	real<lower=0> lP1;
	real<lower=0> lP2; 
	real<upper=0> ER24_20;
	real<lower=0> k600;
	real<lower=0> sigma;
}
transformed parameters {
	vector [maxTime] gpp = rep_vector(0, maxTime);
	vector [maxTime] er = rep_vector(0, maxTime);
	vector [maxTime] DO_pr; 

	DO_pr[1] = DOinitial;
	for(i in 2:maxTime) {
		real ddodt;
		real rf;
		rf = computeRF(temp[i-1], pressure, DO_pr[i-1], k600);
		gpp[i] = computeGPP(PAR[i-1], lP1, lP2);
		er[i] = computeER(temp[i-1], ER24_20);
		ddodt = (gpp[i] + er[i] + rf) / z;
		DO_pr[i] = DO_pr[i-1] + ddodt * dt;
	}
}
model {
	for(i in 1:nDO) {
		DO[i] ~ normal(DO_pr[time[i]], sigma);
	}
	k600 ~ normal((1162 * pow(slope, 0.77) * pow(velocity, 0.85))/(24*60), 0.0001462944 + 0.0012564499 * slope + 0.0124307051 * velocity + 0.0961094198 * slope * velocity);
	lP1 ~ normal(9, 1);
	lP2 ~ normal(9, 1);
	ER24_20 ~ normal(0, 10);
}


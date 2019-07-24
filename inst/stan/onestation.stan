functions {
	real computeRF(real temp, real pressure, real DO, real k600);
	real osat(real temp, real P);
	real kT(real temp, real k600);
	real computeGPP(real PAR, real lP1, real lP2);
	real computeER(real temp, real ER24_20);
	real k_mu(real slope, real velocity);
	real k_sd(real slope, real velocity);

	#include "functions.stan"
}
data {
	int<lower = 1> nDO; // DO number of observations;

	vector<lower=0> [nDO] DO; // DO Observations
	int<lower = 1> time [nDO]; // time of each DO observation, in minutes
	real<lower = 0> DO_initial; // estimate of DO concentration at time = 1

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
	real lP1;
	real lP2; 
	real<upper=0> ER24_20;
	real<lower=0> k600;
	real<lower=0> sigma;
	real<lower=0> DO_i; // initial value DO
}
transformed parameters {
	vector [maxTime] gpp = rep_vector(0, maxTime);
	vector [maxTime] er = rep_vector(0, maxTime);
	vector [maxTime] DO_pr; 

	DO_pr[1] = DO_i;
	for(i in 2:maxTime) {
		real ddodt;
		real rf;
		rf = computeRF(temp[i-1], pressure, DO_pr[i-1], k600);
		gpp[i] = computeGPP(PAR[i-1], lP1, lP2);
		er[i] = computeER(temp[i-1], ER24_20);
		ddodt = ((gpp[i] + er[i] + rf) / z) / (24*60);
		DO_pr[i] = DO_pr[i-1] + ddodt * dt;
	}
}
model {
	for(i in 1:nDO) {
		DO[i] ~ normal(DO_pr[time[i]], sigma);
	}
	k600 ~ normal(k_mu(slope, velocity), k_sd(slope, velocity));
	lP1 ~ normal(9, 3); // considering the log scale, these are VERY weakly informative
	lP2 ~ normal(9, 3); // considering the log scale, these are VERY weakly informative
	ER24_20 ~ normal(0, 10);
	sigma ~ normal(0, 10);
	DO_i ~ normal(DO_initial, sigma);
}


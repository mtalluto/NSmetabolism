functions {
	real computeRF(real temp, real pressure, real DO, real k600);
	real osat(real temp, real P);
	real kT(real temp, real k600);
	real computeGPP(real PAR, real lP1, real lP2);
	real computeER(real temp, real ER24_20);

	#include "include/functions.stan"
}
data {
	int<lower = 1> nDO; // DO number of observations;

	vector<lower=0> [nDO] DO; // DO Observations
	vector [nDO] temp;
	vector <lower = 0> [nDO] light;
	vector <lower = 0> [nDO] pressure;

	real <lower = 0> depth;
	real <lower = 0> delta_t; // length of time (minutes) of a single time step

	// parameters for the priors; mean and scale for each
	vector [2] lp1_pr;
	vector [2] lp2_pr;
	vector [2] k_pr;
	vector [2] gpp_pr;
	vector [2] er_pr;
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
	vector [nDO] DO_pr;
	real gpp = 0;
	real er = 0;

	DO_pr[1] = DO_i;
	{
		vector [nDO] gpp_t = rep_vector(0, nDO);
		vector [nDO] er_t = rep_vector(0, nDO);
		for(i in 2:nDO) {
			real ddodt;
			real rf;
			rf = computeRF(temp[i-1], pressure[i-1], DO_pr[i-1], k600);
			gpp_t[i] = computeGPP(light[i-1], lP1, lP2);
			er_t[i] = computeER(temp[i-1], ER24_20);
			ddodt = ((gpp_t[i] + er_t[i] + rf) / depth) / (24*60);
			DO_pr[i] = DO_pr[i-1] + ddodt * delta_t;
			gpp += delta_t * (gpp_t[i] / (24*60));
			er += delta_t * (er_t[i] / (24*60));
		}
	}
}
model {
	DO ~ normal(DO_pr, sigma);
	k600 ~ normal(k_pr[1], k_pr[2]);
	lP1 ~ normal(lp1_pr[1], lp1_pr[2]);
	lP2 ~ normal(lp2_pr[1], lp2_pr[2]);
	ER24_20 ~ normal(0, 10);
	sigma ~ normal(0, 2); // this is on the scale of DO measurements, so pretty small
	DO_i ~ normal(DO_i, sigma);

	// these raise Jacobian warnings
	// however these are likelihood statements, so the jacobian correction isn't really needed
	// see https://discourse.mc-stan.org/t/non-invertible-change-of-variables/1962/22
	gpp ~ normal(gpp_pr[1], gpp_pr[2]);
	er ~ normal(er_pr[1], er_pr[2]);
}


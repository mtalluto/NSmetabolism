functions {
	real kT(real temp, real k600);
	// real computeRF(real temp, real pressure, real DO, real k600);
	// real osat(real temp, real P);
	// real computeAdvection(real inputDO, real outputDO, real Q, real area, real dx);
	// real computeGPP(real PAR, real lP1, real lP2);
	// real computeER(real temp, real ER24_20);

	#include "functions.stan"
	// real kT(real temp, real k600) {
	// 	real Sc;
	// 	// compute Schmidt number for oxygen
	// 	// parameters from Wanninkhof 1992. appendix
	// 	Sc = 1800.6 - 120.10 * temp + 3.7818 * temp^2 - 0.047608 * temp^3;
		
	// 	// Van de Bogert et al eqn 5
	// 	return k600 * (Sc / 600)^-0.5;
	// }
/*
	// pressure units must be in hPa
	real pressureCorrection (real P, real elev, real newElev) {
		real a = 2.25577e-5;
		real = 5.25588;
		real seaLevelP;

		P *= 100; // convert from hPa
		seaLevelP = P / pow(1 - a * elev, b);
		return((seaLevelP * pow(1 - a * newElev, b))/100);
	}

	// perform inverse distance weighting with a correction for elevation
	// vals: input elevation; assumed to be at 0m elevation
	// dist: distance from focal point to vals
	// elevOut: elevation of focal site
	// n: number of pressure observations
	real idw_pressure(vector vals, vector dist, real elevOut, int n) {
		vector weight = 1/pow(dist, 2);
		for(i in 1:n) {
			vals[i] = pressureCorrection(vals[i], 0, elevOut); // convert back from sea level
		}
		weight = weight / sum(weight);
		return(sum(vals * weight));
	}

	// for all vectors, 1 is downstream, 2:n are upstream
	// vals: measured values
	// nbQ: discharge of neighbor sties
	// dist: distance along river to sites
	// Q: discharge of focal site
	real idw_river (vector [3] vals, vector [3] nbQ, vector[3] dist, real Q) {
		vector[3] QR; // discharge ratio to use; should always have upstream in numerator
		vector[3] weights;
		for(i in 1:3) {
			if(dist[i] == 0)
				return(vals[i])
			QR[i] = nbQ / Q;
		}
		QR[1] = 1/QR[1]; // correct the downstream site to put upstream on numerator

		weights = QR / pow(dist, 2);
		weights = weights / sum(weights);
		return sum(vals * weights);
	}

	// simple linear interpolation between two points
	// x: two x points
	// y: two y points
	// xnew: location of the point to interpolate
	real approx(vector x, vector y, xnew) {
		if(xnew == x[1])
			return(y[1]);
		if(xnew == x[2])
			return(y[2]);
		return(y[1] + (xnew - x[1]) * ((y[2] - y[1])/(x[2] - x[1])));
	}
*/
}

data {
	// sample sizes
	// int<lower = 1> nDO; 	// number of DO observations
	int<lower = 2> maxTime;	// number of time steps
	int<lower = 1> nSites; 	// number of sites where measurements were made

	// DUMMY VARIABLES FOR TESTING
	matrix [nSites, maxTime] dummyKT;


	// site-level variables

	// reach-level variables

	// other-scale variables
	// waterTempMeasured: water temperature is measured with DO, but we interpolate it to every
	// minute and pass in matrix format
	matrix [nSites, maxTime] waterTempMeasured;

	/*
	int<lower = 1> nReaches;
	real<lower=0> dt;		// length (in minutes) of a time step

	// site characteristics
	vector<lower=0> [nSites] DOinitial;
	vector<lower=0> [nSites] Q;
	vector<lower=0> [nSites] area;
	vector<lower=0> [nSites] dx;
	vector<lower=0> [nSites] depth;
	vector [nSites] elevation;
	int<lower=1, upper = nReaches> reachID [nSites];

	// reach characteristics
	vector<lower = 0> [nReaches] slope;
	vector<lower = 0> [nReaches] velocity;

	// measured variables and variables for keeping track of them
	int <lower = 1> nTempSites;
	// for water temperature, for each pixel/site, we keep track of 2 upstream neigbors 
	// (indices 2:3), and a downstream neighbor (index 1); we also track the distance to each
	int <lower = 1, upper = nTempSites> waterTempNbs [nSites, 3] ; 
	matrix <lower = 1, upper = nTempSites> [nSites, 3] waterTempDist;
	// for each site where temperature is measured, we have the measurement, as well
	// as that site's pixelID so we can get back to discharge and other data
	int <lower =1 , upper = nSites> waterTempSiteIDs [nTempSites]; // pointer back to pixelID
	matrix<lower = 0> [nsites, 2] coords;
	vector<lower = 0> [nDO] DO;
	vector<lower = 0> [nDO] DOtimes;
	vector<lower = 0> [nDO] DOsites;
	
	// upstream sites; first column is for the main (larger) upstream pixel
	// second column will only be used for confluences; weight should be 0 for non-confluences
	matrix<lower=0> [nSites, 2] usWeight;
	int<lower = 1, upper = nSites> usNb [nSites, 2]; // neighbor indices for upstream pixels

	// lateral input for each site follows the same structure as upstream sites
	// there is a weight based on discharge and a value
	vector<lower = 0> [nSites] latWeight;
	vector<lower=0> [nSites] latInputDO;

	// data relating to atmospheric pressure
	int<lower = 0> nPressure;
	matrix<lower = 0> [nPressure, 2] prCoords;
	matrix<lower = 0> [nPressure, maxTime] pressure;
	vector [nPressure] prElev;

	// data relating to light
	int <lower = 1> nLightTimes;
	matrix [nSites, nLightTimes] light;
	vector [nLightTimes] lightTimes;

	// adjustable priors
	real logP1_pr_mu; // suggested from 1 station - mean of 9
	real logP1_pr_sd; // suggested from 1 station - sd of 1
	real logP2_pr_mu; 
	real logP2_pr_sd;
*/
}
/*
transformed data {
	matrix <lower = 0> [nPressure, maxTime] pressureSeaLevel;
	matrix <lower = 0> [nSites, nPressure] pressureDist;
	for(i in 1:nSites) {
		for(j in 1:nPressure) {
			pressureDist[i, j] = sqrt(pow(coords[i, 1] - prCoords[j, 1], 2) + 
				pow(coords[i, 2] - prCoords[j, 2], 2))
		}
	}

	for(i in 1:nPressure) {
		for(j in 1:maxTime) {
			pressureSeaLevel[i, j] = pressureCorrection(pressure[i, j], prElev[i], 0)
		}
	}
}
*/
parameters {
	real k600;

	real<lower = 0> sigma;
	/*
	vector<lower=0> [nReaches] k600;
	vector [nReaches] lP1;
	vector [nReaches] lP2;
	vector<lower=0> [nReaches] ER24_20;
	*/
}

transformed parameters {
	matrix [nSites, maxTime] dummyKTTheo;
	for(i in 1:nSites) {
		for(j in 1:maxTime) {
			dummyKTTheo[i,j] = kT(waterTempMeasured[i,j], k600);
		}
	}
	/*
	matrix [nSites, maxTime] gpp = rep_matrix(0, nSites, maxTime);
	matrix [nSites, maxTime] er = rep_matrix(0, nSites, maxTime);
	matrix [nSites, maxTime] DOpr; 
	int ltTimeIndex = 1;

	for(si in 1:nSites) {
		DOpr[si,1] = DOinitial[si];
	}

	for(ti in 2:maxTime) {
		if(lightTimes[ltTimeIndex] > ti)
			ltTimeIndex += 1;
		for(si in 1:nSites) {
			real ddodt;
			real rf;
			real inputDO;
			real adv;
			real waterTemp;
			real pressure;
			real light;
			int reach = reachID[si];
			vector [2] usMeasuredWaterTemp;
			vector [2] usWaterTempQ;
			real dsMeasuredWaterTemp;
			real siPressure;
			real siLight;

			// deal with water temperature interpolation
			// this is a somewhat complicated index lookup so break out the indies a bit
			// for clarity
			{
				int nbIDs [3] = waterTempNbs[si, ];
				int nbPixIDs = waterTempSiteIDs[nbIDs];
				waterTemp = idw_river(waterTempMeasured[nbIDs, ti], Q[nbPixIDs], 
					waterTempDist[si, ], Q[si]);
			}


			// interpolate pressure
			siPressure = idw_pressure(pressureSeaLevel[,ti], pressureDist[si,], 
				elevation[si], nPressure);

			// interpolate light
			siLight = approx(lightTimes[ltTimeIndex:(ltTimeIndex+1)], 
						light[si, ltTimeIndex:(ltTimeIndex+1)], ti)

			// get input DO from upstream pixel(s)
			// note that stan is somewhat inflexible
			// for some sites to have 2 upstream pixels, ALL sites must have 2 upstream pixels
			// this is handled by having a dummy upstream pixel with a weight of 0
			inputDO = usWeight[si,1] * DOpr[usNb[si,1], ti-1] + 
					usWeight[si,2] * DOpr[usNb[si,2], ti-1] + 
					latWeight[si] * latInputDO[si];

			adv = computeAdvection(inputDO, DOpr[si, ti-1], Q[si], area[si], dx[si]);

			// note that all of these components, including ER, are constant at the reach scale
			// finer resolution may be necessary in the future
			rf = computeRF(waterTemp, siPressure, DOpr[si, ti-1], k600[reach]);
			gpp[si, ti] = computeGPP(siLight, lP1[reach], lP2[reach]);
			er[si, ti] = computeER(waterTemp, ER24_20[reach]);
			ddodt = adv + (gpp[si, ti] + er[si, ti] + rf) / depth[si];
			DO_pr[si, ti] = DO_pr[si, ti-1] + ddodt * dt;
		}
	}
	*/
}
model {
	for(i in 1:nSites) {
		for(j in 1:maxTime) {
			// waterTempMeasured[i,j] ~ normal(k600, sigma);
			dummyKT[i,j] ~ normal(dummyKTTheo[i,j], sigma);
		}
	}
	/*
	for(i in 1:nDO) {
		DO[i] ~ normal(DOpr[DOsites[i], DOtimes[i]], sigma);
	}

	for(i in 1:nReaches) {
		k600[i] ~ normal((1162 * pow(slope[i], 0.77) * pow(velocity[i], 0.85))/(24*60),
				0.0001462944 + 0.0012564499 * slope[i] + 0.0124307051 * velocity[i] + 
				0.0961094198 * slope[i] * velocity[i]);
	}
	
	lP1 ~ normal(logP1_pr_mu, logP1_pr_sd);
	lP2 ~ normal(logP2_pr_mu, logP2_pr_sd);
	ER24_20 ~ normal(0, 10);
	*/
}

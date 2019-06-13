functions {
	real kT(real temp, real k600);
	real computeRF(real temp, real pressure, real DO, real k600);
	real osat(real temp, real P);
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

	/**
	 * Compute pressure at an elevation from pressure at a different elevation
	 *
	 * @param P pressure in hPa
	 * @param elev elevation (meters) where pressure was measured
	 * @param newElev elevation (m) at which to compute pressure
	 * @return pressure (hPa) at specified elevation
	*/
	real pressureCorrection (real P, real elev, real newElev) {
		real a = 2.25577e-5;
		real b = 5.25588;
		real seaLevelP;
		real PPa; // pressure in pascals

		PPa = P * 100; // convert from hPa
		seaLevelP = PPa / pow(1 - a * elev, b);
		return((seaLevelP * pow(1 - a * newElev, b))/100);
	}


	// perform inverse distance weighting with a correction for elevation
	// vals: input elevation; assumed to be at 0m elevation
	// dist: distance from focal point to vals
	// elevOut: elevation of focal site
	// n: number of pressure observations
	// real idw_pressure(vector vals, vector dist, real elevOut, int n) {
	// 	vector weight = 1/pow(dist, 2);
	// 	for(i in 1:n) {
	// 		vals[i] = pressureCorrection(vals[i], 0, elevOut); // convert back from sea level
	// 	}
	// 	weight = weight / sum(weight);
	// 	return(sum(vals * weight));
	// }

	/**
	 * Compute discharge and inverse square distance weighted interpolation based on river distance
	 * all vectors are of length 3 with index 1 being downstream and 2:3 upstream
	 *
	 * @param vals measured values
	 * @param nbQ discharge of neighbor sties
	 * @param dist distance (m) along river to sites
	 * @param Q discharge of focal site
	 * @return interpolated value of vals
	*/
	real idw_river (vector vals, vector nbQ, vector dist, real Q) {
		vector[3] QR; // discharge ratio to use; should always have upstream in numerator
		vector[3] weights;
		vector[3] res;

		// check if site is downstream of itself
		if(dist[1] == 0)
			return(vals[1]);
		
		QR = nbQ / Q;
		QR[1] = 1/QR[1]; // correct the downstream site to put upstream on numerator

		for(i in 1:3)
			weights[i] = QR[i] / pow(dist[i], 2);

		weights = weights / sum(weights);
		for(i in 1:3)
			res[i] = vals[i] * weights[i];
		return sum(res);
	}

	// simple linear interpolation between two points
	// x: two x points
	// y: two y points
	// xnew: location of the point to interpolate
	// real approx(vector x, vector y, xnew) {
	// 	if(xnew == x[1])
	// 		return(y[1]);
	// 	if(xnew == x[2])
	// 		return(y[2]);
	// 	return(y[1] + (xnew - x[1]) * ((y[2] - y[1])/(x[2] - x[1])));
	// }

}

data {
	// sample sizes
	int<lower = 1> nDO; 	// number of DO observations
	int<lower = 2> maxTime;	// number of time steps
	int<lower = 1> nSites; 	// number of sites where measurements were made
	int<lower = 0> nPressure;
	int<lower = 1> nPixels; // number of pixels in the watershed
	int<lower = 1> nReaches;

	// site-level variables

	// reach-level variables



	// real<lower=0> dt;		// length (in minutes) of a time step

	// // site characteristics
	// vector<lower=0> [nSites] DOinitial;
	vector<lower=0> [nPixels] Q;
	// vector<lower=0> [nSites] area;
	// vector<lower=0> [nSites] dx;
	// vector<lower=0> [nSites] depth;
	// vector [nSites] elevation;
	// int<lower=1, upper = nReaches> reachID [nSites];

	// // reach characteristics
	// vector<lower = 0> [nReaches] slope;
	// vector<lower = 0> [nReaches] velocity;

	/*
		WATER TEMPERATURE
	*/
	// waterTempMeasured: water temperature is measured with DO, but we interpolate it to every
	// minute and pass in matrix format
	matrix [nSites, maxTime] waterTempMeasured;
	// for water temperature, for each pixel/site, we keep track of 2 upstream neigbors 
	// (indices 2:3), and a downstream neighbor (index 1); we also track the distance to each
	int <lower = 1, upper = nSites> waterTempNbs [nPixels, 3] ; 
	matrix <lower = 1, upper = nSites> [nPixels, 3] waterTempDist;
	// for each site where temperature is measured, we have the measurement, as well
	// as that site's pixelID so we can get back to discharge and other data
	int <lower =1 , upper = nPixels> waterTempSiteIDs [nSites]; // pointer back to pixelID



	// matrix<lower = 0> [nsites, 2] coords;
	vector<lower = 0> [nDO] DO;
	int<lower = 0, upper = maxTime> DOtimes [nDO];
	int<lower = 1, upper = nPixels> DOpixels [nDO];
	
	// // upstream sites; first column is for the main (larger) upstream pixel
	// // second column will only be used for confluences; weight should be 0 for non-confluences
	// matrix<lower=0> [nSites, 2] usWeight;
	// int<lower = 1, upper = nSites> usNb [nSites, 2]; // neighbor indices for upstream pixels

	// // lateral input for each site follows the same structure as upstream sites
	// // there is a weight based on discharge and a value
	// vector<lower = 0> [nSites] latWeight;
	// vector<lower=0> [nSites] latInputDO;



	// data relating to atmospheric pressure
	//matrix<lower = 0> [nPressure, 2] prCoords;
	matrix<lower = 0> [nPressure, maxTime] pressure; // atmospheric pressure, in hPa (i.e., mbar)
	vector [nPressure] prElev;

/*
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

transformed data {
	matrix <lower = 0> [nPressure, maxTime] pressureSeaLevel;
	// matrix <lower = 0> [nSites, nPressure] pressureDist;
	// for(i in 1:nSites) {
	// 	for(j in 1:nPressure) {
	// 		pressureDist[i, j] = sqrt(pow(coords[i, 1] - prCoords[j, 1], 2) + 
	// 			pow(coords[i, 2] - prCoords[j, 2], 2))
	// 	}
	// }

	for(i in 1:nPressure) {
		for(j in 1:maxTime) {
			pressureSeaLevel[i, j] = pressureCorrection(pressure[i, j], prElev[i], 0);
		}
	}
}

parameters {

	real<lower = 0> sigma;
	vector<lower=0> [nReaches] k600;
	// vector [nReaches] lP1;
	// vector [nReaches] lP2;
	// vector<lower=0> [nReaches] ER24_20;

}

transformed parameters {	
	// matrix [nPixels, maxTime] gpp = rep_matrix(0, nPixels, maxTime);
	// matrix [nPixels, maxTime] er = rep_matrix(0, nPixels, maxTime);
	matrix [nPixels, maxTime] DOpr; 
	// int ltTimeIndex = 1;

	// for(si in 1:nSites) {
	// 	DOpr[si,1] = DOinitial[si];
	// }

	for(ti in 2:maxTime) {
		// if(lightTimes[ltTimeIndex] > ti)
			// ltTimeIndex += 1;
		for(pix in 1:nPixels) {
			// real ddodt;
			real rf;
			// real inputDO;
			// real adv;
			real waterTemp;
			// real pressure;
			// real light;
			int reach;
			// int reach = reachID[si];
			// vector [2] usMeasuredWaterTemp;
			// vector [2] usWaterTempQ;
			// real dsMeasuredWaterTemp;
			real pixPressure;
			// real siLight;

			/* 
				deal with water temperature interpolation
				this is a somewhat complicated index lookup so break out the indies a bit
				for clarity
			*/
			{
				int nbIDs [3] = waterTempNbs[pix, ];
				int nbPixIDs [3] = waterTempSiteIDs[nbIDs];
				waterTemp = idw_river(waterTempMeasured[nbIDs, ti], Q[nbPixIDs], 
					to_vector(waterTempDist[pix, ]), Q[pix]);
			}


			// interpolate pressure
			// siPressure = idw_pressure(pressureSeaLevel[,ti], pressureDist[si,], 
			// 	elevation[si], nPressure);

			// interpolate light
			// siLight = approx(lightTimes[ltTimeIndex:(ltTimeIndex+1)], 
			// 			light[si, ltTimeIndex:(ltTimeIndex+1)], ti)

			// get input DO from upstream pixel(s)
			// note that stan is somewhat inflexible
			// for some sites to have 2 upstream pixels, ALL sites must have 2 upstream pixels
			// this is handled by having a dummy upstream pixel with a weight of 0
			// inputDO = usWeight[si,1] * DOpr[usNb[si,1], ti-1] + 
			// 		usWeight[si,2] * DOpr[usNb[si,2], ti-1] + 
			// 		latWeight[si] * latInputDO[si];

			// adv = computeAdvection(inputDO, DOpr[si, ti-1], Q[si], area[si], dx[si]);

			// note that all of these components, including ER, are constant at the reach scale
			// finer resolution may be necessary in the future

			// BEGIN GARBAGE TESTING CODE
				pixPressure = pressureSeaLevel[1,ti];
				DOpr[pix, ti-1] = DO[1];
				reach = 1;
			// END GARBAGE

			rf = computeRF(waterTemp, pixPressure, DOpr[pix, ti-1], k600[reach]);
			// gpp[si, ti] = computeGPP(siLight, lP1[reach], lP2[reach]);
			// er[si, ti] = computeER(waterTemp, ER24_20[reach]);
			// ddodt = adv + (gpp[si, ti] + er[si, ti] + rf) / depth[si];
			// DO_pr[si, ti] = DO_pr[si, ti-1] + ddodt * dt;
		}
	}
}
model {

	// for(i in 1:nDO) {
	// 	DO[i] ~ normal(DOpr[DOpixels[i], DOtimes[i]], sigma);
	// }
/*
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

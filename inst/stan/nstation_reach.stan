functions {
	real kT(real temp, real k600);
	real computeRF(real temp, real pressure, real DO, real k600);
	real osat(real temp, real P);
	real computeAdvection(real inputDO, real outputDO, real Q, real area, real dx);
	real computeGPP(real PAR, real lP1, real lP2);
	real computeER(real temp, real ER24_20);

	#include "functions.stan"


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

	/**
	 * compute a weighted average
	 *
	 * @param vals Values to average
	 * @param weights Weights to assign
	 * @param n Length of vectors
	*/
	real weighted_avg (vector vals, vector weights, int n) {
		vector [n] normWts;
		vector [n] res;
		normWts = weights / sum(weights);
		for(i in 1:n)
			res[i] = vals[i] * normWts[i];
		return sum(res);
	}


	/**
	 * perform inverse distance weighting with a correction for elevation
	 *
	 * @param vals input elevation; assumed to be at 0m elevation
	 * @param dist distance from focal point to vals
	 * @param elevOut elevation of focal site
	 * @param n number of pressure observations
	 * @return interpolated value of vals
	*/
	real idw_pressure(vector vals, vector dist, real elevOut, int n) {
		vector [n] weight;
		vector [n] valCorr;
		vector [n] res;

		for(i in 1:n) {
			weight[i] = 1/pow(dist[i], 2);
			valCorr[i] = pressureCorrection(vals[i], 0, elevOut); // convert back from sea level
		}
		return weighted_avg(vals, weight, n);
	}



	/**
	 * Compute discharge and inverse square distance weighted interpolation based on river distance
	 * all vectors are of length 3 with index 1 being downstream and 2:3 upstream
	 *
	 * @param vals measured values
	 * @param nbQ discharge of neighbor sties
	 * @param dist distance (m) along river to sites
	 * @param Q discharge of focal site
	 * @return 
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

		return weighted_avg(vals, weights, 3);
	}

	/**
	 * simple linear interpolation between two points
	 *
	 * @param x two x points
	 * @param y two y points
	 * @param xnew location of the point to interpolate
	*/
	real approx(vector x, vector y, real xnew) {
		if(xnew == x[1])
			return(y[1]);
		if(xnew == x[2])
			return(y[2]);
		return(y[1] + (xnew - x[1]) * ((y[2] - y[1])/(x[2] - x[1])));
	}

}

data {
	// sample sizes
	int<lower = 1> nDO; 	// number of DO observations
	int<lower = 2> maxTime;	// number of time steps
	int<lower = 1> nSites; 	// number of sites where measurements were made
	int<lower = 0> nPressure;
	int<lower = 1> nPixels; // number of pixels in the watershed
	int<lower = 1> nReaches;

	real<lower=0> dt;		// length (in minutes) of a time step

	// pixel characteristics
	vector<lower=0> [nPixels] DOinitial;
	vector<lower=0> [nPixels] Q;
	vector<lower=0> [nPixels] area;
	vector<lower=0> [nPixels] dx;
	vector<lower=0> [nPixels] depth;
	vector [nPixels] elevation;
	int<lower=1, upper = nReaches> reachID [nPixels];
	matrix [nPixels, 2] coords;

	// reach characteristics
	vector<lower = 0> [nReaches] slope;
	vector<lower = 0> [nReaches] velocity;

	/*
		WATER TEMPERATURE
	*/
	// waterTempMeasured: water temperature is measured with DO, but we interpolate it to every
	// minute and pass in matrix format
	matrix [nSites, maxTime] waterTempMeasured;
	// for water temperature, for each pixel/site, we keep track of 2 upstream neigbors 
	// (indices 2:3), and a downstream neighbor (index 1); we also track the distance to each
	int <lower = 1, upper = nPixels> waterTempNbs [nPixels, 3] ; 
	matrix [nPixels, 3] waterTempDist;
	// for each site where temperature is measured, we have the measurement, as well
	// as that site's pixelID so we can get back to discharge and other data
	int <lower =1 , upper = nPixels> waterTempSiteIDs [nSites]; // pointer back to pixelID
	int <lower = 1, upper = nSites> waterTempIndices [nPixels, 3] ;


	vector<lower = 0> [nDO] DO;
	int<lower = 0, upper = maxTime> DOtimes [nDO];
	int<lower = 1, upper = nPixels> DOpixels [nDO];
	
	// upstream sites; first column is for the main (larger) upstream pixel
	// second column will only be used for confluences; weight should be 0 for non-confluences
	matrix<lower=0> [nPixels, 2] usWeight;
	int<lower = 1, upper = nPixels> usNb [nPixels, 2]; // neighbor indices for upstream pixels

	// lateral input for each site follows the same structure as upstream sites
	// there is a weight based on discharge and a value
	vector [nPixels] latWeight;
	vector<lower=0> [nPixels] latInputDO;



	/*
		PRESSURE
	*/
	matrix [nPressure, 2] prCoords;
	matrix<lower = 0> [nPressure, maxTime] pressure; // atmospheric pressure, in hPa (i.e., mbar)
	vector [nPressure] prElev;

	/*
		LIGHT
	*/
	int <lower = 1> nLightTimes;
	matrix [nPixels, nLightTimes] light;
	vector [nLightTimes] lightTimes;

	// adjustable priors
	real logP1_pr_mu; // suggested from 1 station - mean of 9
	real logP1_pr_sd; // suggested from 1 station - sd of 1
	real logP2_pr_mu; 
	real logP2_pr_sd;
}

transformed data {
	matrix <lower = 0> [nPressure, maxTime] pressureSeaLevel;
	matrix <lower = 0> [nPixels, nPressure] pressureDist;
	for(i in 1:nPixels) {
		for(j in 1:nPressure) {
			pressureDist[i, j] = sqrt(pow(coords[i, 1] - prCoords[j, 1], 2) + 
				pow(coords[i, 2] - prCoords[j, 2], 2));
		}
	}

	for(i in 1:nPressure) {
		for(j in 1:maxTime) {
			pressureSeaLevel[i, j] = pressureCorrection(pressure[i, j], prElev[i], 0);
		}
	}
}

parameters {

	real<lower = 0> sigma;
	vector<lower=0> [nReaches] k600;
	vector [nReaches] lP1;
	vector [nReaches] lP2;
	vector<lower=0> [nReaches] ER24_20;

}

transformed parameters {	

}
model {
	{
	matrix [nPixels, maxTime] gpp = rep_matrix(0, nPixels, maxTime);
	matrix [nPixels, maxTime] er = rep_matrix(0, nPixels, maxTime);
	matrix [nPixels, maxTime] DOpr; 

	for(pix in 1:nPixels) {
		DOpr[pix,1] = DOinitial[pix];
	}
	{
		int ltTimeIndex = 1;

		for(ti in 2:maxTime) {
			real ctime = ti;
			if(lightTimes[ltTimeIndex] > ti)
				ltTimeIndex += 1;
			for(pix in 1:nPixels) {
				real ddodt;
				real rf;
				real inputDO;
				real adv;
				real waterTemp;
				int reach = reachID[pix];
				real pixPressure;
				real pixLight;

				/* 
					deal with water temperature interpolation
					this is a somewhat complicated index lookup so break out the indies a bit
					for clarity
				*/
				{
					int nbIDs [3] = waterTempNbs[pix, ];
					waterTemp = idw_river(waterTempMeasured[waterTempIndices[pix,], ti], 
						Q[nbIDs], to_vector(waterTempDist[pix, ]), Q[pix]);
				}


				// interpolate pressure
				pixPressure = idw_pressure(pressureSeaLevel[,ti], to_vector(pressureDist[pix,]), 
					elevation[pix], nPressure);

				// interpolate light
				pixLight = approx(lightTimes[ltTimeIndex:(ltTimeIndex+1)], 
							to_vector(light[pix, ltTimeIndex:(ltTimeIndex+1)]), ctime);



				// get input DO from upstream pixel(s)
				// note that stan is somewhat inflexible
				// for some sites to have 2 upstream pixels, ALL sites must have 2 upstream pixels
				// this is handled by having a dummy upstream pixel with a weight of 0
				inputDO = usWeight[pix,1] * DOpr[usNb[pix,1], ti-1] + 
						usWeight[pix,2] * DOpr[usNb[pix,2], ti-1] + 
				 		latWeight[pix] * latInputDO[pix];
				adv = computeAdvection(inputDO, DOpr[pix, ti-1], Q[pix], area[pix], dx[pix]);
				print(adv);
				// note that all of these components, including ER, are constant at the reach scale
				// finer resolution may be necessary in the future
				rf = computeRF(waterTemp, pixPressure, DOpr[pix, ti-1], k600[reach]);
				print(rf);
				gpp[pix, ti] = computeGPP(pixLight, lP1[reach], lP2[reach]);
				print(gpp[pix, ti]);
				er[pix, ti] = computeER(waterTemp, ER24_20[reach]);
				print(er[pix, ti]);
				ddodt = adv + (gpp[pix, ti] + er[pix, ti] + rf) / depth[pix];
				DOpr[pix, ti] = DOpr[pix, ti-1] + ddodt * dt;
			}
		}
	}
	for(i in 1:nDO) {
		DO[i] ~ normal(DOpr[DOpixels[i], DOtimes[i]], sigma);
	}
}
	for(i in 1:nReaches) {
		k600[i] ~ normal((1162 * pow(slope[i], 0.77) * pow(velocity[i], 0.85))/(24*60),
				0.0001462944 + 0.0012564499 * slope[i] + 0.0124307051 * velocity[i] + 
				0.0961094198 * slope[i] * velocity[i]);
	}
	
	lP1 ~ normal(logP1_pr_mu, logP1_pr_sd);
	lP2 ~ normal(logP2_pr_mu, logP2_pr_sd);
	ER24_20 ~ normal(0, 10);

}

functions {
	real computeRF() {

	}

	real computeGPP() {

	}

	real computeER() {

	}

	real idw() {

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

	// performs nearest neighbor linear interpolation
	real approx(real x1, real y1, real x2, real y2, real xnew) {

	}

	real computeAdvection(real inputDO, real outputDO, real Q, real area, real dx) {
		real inputMass = Q * inputDO;
		real outputMass = Q * outputDO;
		return (-1/area) * (outputMass - inputMass)/dx;
	}
}

data{
	int<lower=1> nSites; 	// number of sites
	int<lower=2> maxTime;	// number of time steps
	int<lower = 1> nReaches;
	real<lower=0> dt;		// length (in minutes) of a time step

	// site characteristics
	vector<lower=0> [nSites] DOinitial;
	vector<lower=0> [nSites] Q;
	vector<lower=0> [nSites] area;
	vector<lower=0> [nSites] dx;
	vector<lower=0> [nSites] depth;
	int<lower=1, upper = nReaches> reachID [nSites];

	// measured variables and variables for keeping track of them
	int <lower = 1> nTempSites;
	// for water temperature, for each pixel/site, we keep track of 2 upstream neigbors 
	// (indices 2:3), and a downstream neighbor (index 1); we also track the distance to each
	int <lower = 1, upper = nTempSites> waterTempNbs [nSites, 3] ; 
	matrix <lower = 1, upper = nTempSites> [nSites, 3] waterTempDist;
	// for each site where temperature is measured, we have the measurement, as well
	// as that site's pixelID so we can get back to discharge and other data
	int <lower =1 , upper = nSites> waterTempSiteIDs [nTempSites]; // pointer back to pixelID
	matrix [nTempSites, maxTime] waterTempMeasured;
	
	// upstream sites; first column is for the main (larger) upstream pixel
	// second column will only be used for confluences; weight should be 0 for non-confluences
	matrix<lower=0> [nSites, 2] usWeight;
	int<lower = 1, upper = nSites> usNb [nSites, 2]; // neighbor indices for upstream pixels

	// lateral input for each site follows the same structure as upstream sites
	// there is a weight based on discharge and a value
	vector<lower = 0> [nSites] latWeight;
	vector<lower=0> [nSites] latInputDO;
}
parameters {
	vector<lower=0> [nReaches] k600;
	vector<lower=0> [nReaches] lP1;
	vector<lower=0> [nReaches] lP2;
	vector<lower=0> [nReaches] ER24_20;
}

transformed parameters {
	matrix [nSites, maxTime] gpp = rep_matrix(0, nSites, maxTime);
	matrix [nSites, maxTime] er = rep_matrix(0, nSites, maxTime);
	matrix [nSites, maxTime] DOpr; 

	for(si in 1:nSites) {
		DOpr[si,1] = DOinitial[si];
	}

	for(ti in 2:maxTime) {
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

			// deal with water temperature interpolation
			// this is a somewhat complicated index lookup so break out the indies a bit
			// for clarity
			{
				int nbIDs [3] = waterTempNbs[si, ];
				int nbPixIDs = waterTempSiteIDs[nbIDs];
				waterTemp = idw_river(waterTempMeasured[nbIDs, ti], Q[nbPixIDs], 
					waterTempDist[si, ], Q[si]);
			}

			// deal with pressure interpolation; same deal as water temp; the indices get complex
			{
				// pressure = idw(...);
			}

			// deal with light interpolation; same deal as water temp; the indices get complex
			{
				// light = approx(...);
			}

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
			rf = computeRF(waterTemp, pressure, DOpr[si, ti-1], k600[reach]);
			gpp[si, ti] = computeGPP(light, lP1[reach], lP2[reach]);
			er[si, ti] = computeER(waterTemp, ER24_20[reach]);
			ddodt = adv + (gpp[si, ti] + er[si, ti] + rf) / depth[si];
			DO_pr[si, ti] = DO_pr[si, ti-1] + ddodt * dt;
		}
	}
}

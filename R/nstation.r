library("WatershedTools")

"	DATA TO DEAL WITH
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
"

"	functions to test
	real kT(real temp, real k600);
	real computeRF(real temp, real pressure, real DO, real k600);
	real osat(real temp, real P);
	real computeAdvection(real inputDO, real outputDO, real Q, real area, real dx);
	real computeGPP(real PAR, real lP1, real lP2);
	real computeER(real temp, real ER24_20);

	real pressureCorrection (real P, real elev, real newElev) {
	real idw_pressure(vector vals, vector dist, real elevOut, int n) {
	real idw_river (vector [3] vals, vector [3] nbQ, vector[3] dist, real Q) {
	real approx(vector x, vector y, xnew) {


	ALSO TEST 
		transformed_data (can do in line with functions where appropriate)
		transformed_parameters
		finally model
"
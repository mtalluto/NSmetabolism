// stan file included for historical purposes; do not calibrate with stan (too slow)

functions {
	real osat(real temp, real BP) {
		real DOo;
		real theta;
		real mu;
		real Fp;
		
		DOo = exp(-139.34411+(157570.1/(temp + 273.15)) - (66423080/(temp + 273.15)^2) + (12438000000/(temp + 273.15)^3) - (862194900000/(temp + 273.15)^4));
		theta = 0.000975 - 0.00001426*temp + 0.00000006436*temp^2;
		mu = exp(11.8571 - (3840.7/(temp + 273.15)) - (216961/(temp + 273.15)^2));
		Fp = ((BP - mu)*(1-theta*BP))/((1-mu)*(1-theta));
		return DOo*Fp;
	}

	real kCor (real temp, real K600) {
		real Sc;

		Sc = 1568 - 86.04*temp + 2.142*temp^2 - 0.0216*temp^3;
		return K600*(Sc/600)^-0.5;
	}
}
data {

	int<lower = 1> nt; // number of time steps;
	real<lower=0> dt; // size of time step
	vector<lower = 0> [nt] PAR;
	vector [nt] temp;
	real do_initial;
	real<lower=0> z;
	real BP;
	vector<lower=0> [nt] DO;
}
parameters {
	real<lower=0> P1;
	real<lower=0> P2; 
	real<upper=0> ER24;
	real<lower=0> K600;
	real<lower=0> sigma;
}
transformed parameters {
	vector [nt-1] production;
	vector [nt-1] er;
	vector [nt] do_pr;
	do_pr[1] = do_initial;
	for(t in 2:nt)
	{
		real reaeration;
		real ddo;
		production[t-1] = PAR[t-1] / (P1 + P2*PAR[t-1]);
		er[t-1] = (ER24 / (60*24)) * (1.045 ^ (temp[t-1] - 20));
		reaeration = kCor(temp[t-1], K600) * (osat(temp[t-1], BP) - DO[t-1]);
		ddo = (production[t-1] + er[t-1] + reaeration) / z;
		do_pr[t] = do_pr[t-1] + ddo*dt;
	}
}
model {
	DO ~ normal(do_pr, sigma);
}
generated quantities {
	real gpp24;
	real insituER24;
	gpp24 = sum(production);
	insituER24 = sum(er);
}

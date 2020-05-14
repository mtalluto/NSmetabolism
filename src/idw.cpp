#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
using std::vector;

//' Inverse distance weighted interpolation
//'
//' @details Interpolation is done assuming locations are in the rows of vals. Thus, if vals
//' is an n by m matrix, each of the m output values will be the inverse distance weighted 
//' average for the n sites
//'
//' @param vals A location-by-value matrix
//' @param dist Distance to each value
//' @param pow Power to use for distance
//' @return The inverse distance weighted average of vals
//' @export
// [[Rcpp::export]]
NumericVector idw_matrix(const NumericMatrix &vals, NumericVector dist, double pow = 2) {
	if(vals.nrow() != dist.size())
		stop("length(dist) must equal nrow(vals)");

	// convert distances to weights
	std::vector<double> wt;
	{
		double sumwt = 0;
		for(int i = 0; i < dist.size(); ++i) {
			wt.push_back(1/std::pow(dist[i], pow));
			sumwt += wt[i];
		}
		for(auto &w : wt)
			w = w / sumwt;
	}
	
	NumericVector result (vals.ncol(), 0.0);

	for(int i = 0; i < vals.ncol(); ++i) {
		for(int j = 0; j < vals.nrow(); ++j)
			result[i] += wt[j] * vals(j,i);
	}
	return result;
}

/**
	*' Compute the average of x, weighted by w
*/
double weighted_avg(const vector<double> &x, const vector<double> &w) {
	double sum = std::accumulate(w.cbegin(), w.cend(), 0.0);
	vector<double> res;

	for(unsigned i = 0; i < w.size(); ++i) {
		double normW = w.at(i) / sum;
		res.push_back(x.at(i) * normW);
	}

	return std::accumulate(res.begin(), res.end(), 0.);
}


/**
 * @brief      Perform inverse distance weighting with a 
 * 				correction for discharge for a single point
 *
 * @param[in]  vals  input values of variable to interpolate
 * @param[in]  dist  distance from focal point to vals; negative for downstream sites
 * @param[in]  nbQ   discharge of measured sites
 * @param[in]  Q     discharge of focal site
 *
 * @return     interpolated value of vals
 */
double idw_river_pt(const vector<double> &vals, const vector<double> &dist, 
		const vector<double> &nbQ, double Q) {
	// check if the focal site is itentical to any of the sites
	// vector<double>::const_iterator it = std::find(dist.begin(), dist.end(), 0.);
	// if(it != dist.end())
	// 	return(vals.at(it - dist.begin()));

	// compute discharge ratios, with downstream sites always in the denominator
	vector<double> weights;
	for(int i = 0; i < nbQ.size(); ++i) {
		if(dist.at(i) == 0)
			return(vals.at(i));
		double QR = dist.at(i) < 0 ? Q/nbQ.at(i) : nbQ.at(i) / Q;
		weights.push_back(QR / std::pow(dist.at(i), 2));
	}
	return weighted_avg(vals, weights);
}


//' Perform inverse distance weighting with a correction for discharge
//'
//' This function works on a vector of k output sites, where the length
//' of each input list and the input vector Q must be k.
//' For each input site, the corresponding input lists contain vectors of neighbor sites;
//' so for site k with discharge Q[k], vals[[k]] is a vector for the n neighbors of site k.
//'
//' @param vals input values of variable to interpolate
//' @param dist distance from focal point to vals; negative for downstream sites
//' @param nbQ discharge of measured sites
//' @param Q discharge of focal site
//' @return interpolated value of vals
//' @export
// [[Rcpp::export]]
NumericVector idw_river(List vals, List dist, List nbQ, NumericVector Q) {
	if(!(vals.length() == dist.length() & vals.length() == nbQ.length() & 
			vals.length()== Q.length()))
		throw std::length_error("Argument lengths must all match");
	NumericVector result;
	for(int i = 0; i < vals.length(); ++i)
		result.push_back(idw_river_pt(vals(i), dist(i), nbQ(i), Q(i)));
	return result;
}

#include <vector>
#include <Rcpp.h>

void check_Q(const Rcpp::NumericVector &Q);
void check_site_data(const Rcpp::NumericVector &site);
void check_positive(const Rcpp::NumericVector &dat, std::string par);

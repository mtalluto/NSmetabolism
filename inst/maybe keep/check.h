#include <vector>
#include <Rcpp.h>

void check_DO(const Rcpp::NumericVector &DO);
void check_Q(const Rcpp::NumericVector &Q);
void check_pressure(const Rcpp::NumericVector &P);
void check_site_data(const Rcpp::NumericVector &site);
void check_k(const Rcpp::NumericVector &k);
void check_positive(const Rcpp::NumericVector &dat, std::string par);

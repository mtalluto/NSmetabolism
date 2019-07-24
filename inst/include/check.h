#ifndef CHECK_H
#define CHECK_H

#include <vector>
#include <Rcpp.h>
#include <string>

void check_DO(const Rcpp::NumericVector &DO);
void check_Q(const Rcpp::NumericVector &Q);
void check_pressure(const Rcpp::NumericVector &P);
void check_site_data(const Rcpp::NumericVector &site);
void check_light(const Rcpp::NumericVector &light);
void check_k(const Rcpp::NumericVector &k);
void check_ge_zero(const Rcpp::NumericVector &dat, std::string par);
void check_positive(const Rcpp::NumericVector &dat, std::string par);

#endif
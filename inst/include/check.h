#ifndef check_h
#define check_h

#include <Rcpp.h>
#include <string>

void check_light(const Rcpp::NumericVector &light);
void check_ge_zero(const Rcpp::NumericVector &dat, std::string par);
void check_DO(const Rcpp::NumericVector &DO);
void check_pressure(const Rcpp::NumericVector &P);
void check_k(const Rcpp::NumericVector &k);


#endif

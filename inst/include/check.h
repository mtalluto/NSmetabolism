#ifndef check_h
#define check_h

#include <Rcpp.h>
#include <string>

void check_light(const Rcpp::NumericVector &light);
void check_ge_zero(const Rcpp::NumericVector &dat, std::string par);


#endif
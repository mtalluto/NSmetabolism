#include <algorithm>
#include <cmath>
#include <Rcpp.h>
#include "../inst/include/funcs.h"

long double log_sum_exp(long double v1, long double v2) {
	long double mv = std::max(v1, v2);
	return 	mv + std::log(std::exp(v1 - mv) + std::exp(v2 - mv));
}

#ifndef FUNCS_H
#define FUNCS_H

/**
 * @brief      Returns the log of the sum of exponentials
 *
 * @param	v1    First value to sum
 * @param	v2    Second value to sum
 *
 * @return     { the log of the sum of the exponential of v1 and v2: log(exp(v1) + exp(v2)) }
 */
long double log_sum_exp(long double v1, long double v2);


/**
 * Compute pressure at an elevation from pressure at a different elevation
 *
 * @param P pressure in hPa
 * @param elev elevation (meters) where pressure was measured
 * @param newElev elevation (m) at which to compute pressure
 * @return pressure (hPa) at specified elevation
*/
double pressureCorrection (double P, double elev, double newElev);

#endif

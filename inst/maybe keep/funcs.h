#ifndef FUNCS_H
#define FUNCS_H



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

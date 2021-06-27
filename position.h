/** @file position.h
 *  @brief Specification of position
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef POSITION_h_
#define POSITION_h_

/** Position vector (r [m]) from geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
 *  @param [in] lon longitude
 *  @param [in] lat latitude
 *  @param [in] h altitude
 *  @param [out] double
 */
double *position(double lon, double lat, double h);

#endif

/** @file Geodetic.h
 *  @brief Specification of Geodetic
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef GEODETIC_H
#define GEODETIC_H

/** geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]) from given position vector (r [m])
 *  @param [in] r
 *  @param [out] lon
 *  @param [out] lat
 *  @param [out] h
 */
void Geodetic(double *r, double *lon, double *lat, double *h);

#endif

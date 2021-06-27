/** @file AzElPa.h
 *  @brief Specification of AzElPa
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef AZELPA_H
#define AZELPA_H

/** Computes azimuth, elevation and partials from local tangent coordinates
 *  @param [in] s      Topocentric local tangent coordinates (East-North-Zenith frame)
 *  @param [out] A      Azimuth [rad]
 *  @param [out] E      Elevation [rad]
 *  @param [out] dAds   Partials of azimuth w.r.t. s
 *  @param [out] dEds   Partials of elevation w.r.t. s
 */
void AzElPa(double *s, double Az, double El, double *dAds, double *dEds);

#endif

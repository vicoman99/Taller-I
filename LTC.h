/** @file LTC.h
 *  @brief Specification of LTC
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef LTC_H
#define LTC_H

/** Transformation from Greenwich meridian system to local tangent coordinates
 *  @param [in] lon      -Geodetic East longitude [rad]
 *  @param [in] lat      -Geodetic latitude [rad]
 *  @return M        -Rotation matrix from the Earth equator and Greenwich meridian to the local tangent (East-North-Zenith) coordinate system
 */
double **LTC(double lon, double lat);

#endif

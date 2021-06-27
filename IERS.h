/** @file IERS.h
 *  @brief Specification of IERS
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef IERS_H
#define IERS_H

/** Management of IERS time and polar motion data
 *  @param [in] Mjd_UTC
 *  @param [in] interp
 *  @param [out] x_pole
 *  @param [out] y_pole
 *  @param [out] UT1_UTC
 *  @param [out] LOD
 *  @param [out] dpsi
 *  @param [out] deps
 *  @param [out] dx_pole
 *  @param [out] dy_pole
 *  @param [out] TAI_UTC
 */
void IERS(double Mjd_UTC, char interp, double *x_pole, double *y_pole, double *UT1_UTC, double *LOD, double *dpsi, double *deps, double *dx_pole, double *dy_pole, double *TAI_UTC);

#endif/* IERS_H */
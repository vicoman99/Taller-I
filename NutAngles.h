/** @file NutAngles.h
 *  @brief Specification of NutAngles
 *
 *  @author Victor Coman
 *  @date May 2021
 */

#ifndef NUTANGELS_H
#define NUTANGELS_H

/** Nutation in longitude and obliquity
 *  @param [in] Mjd_TT     Modified Julian Date (Terrestrial Time)
 *  @param [out] dpsi,deps  Nutation Angles
 */
void NutAngles(double Mjd_TT, double *dpsi, double *deps);

#endif
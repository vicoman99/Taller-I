/** @file gast.h
 *  @brief Specification of gast
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef GAST_H
#define GAST_H

/** Greenwich Apparent Sidereal Time
 *  @param [in] Mjd_UT1   Modified Julian Date UT1
 *  @return gstime    GAST in [rad]
 */
double gast(double Mjd_UT1);

#endif
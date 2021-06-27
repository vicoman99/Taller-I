/** @file JPL_Eph_DE430.h
 *  @brief Specification of JPL_Eph_DE430
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef JPL_EPH_DE430_H
#define JPL_EPH_DE430_H

/** Computes the sun, moon, and nine major planets equatorial position using JPL Ephemerides
 *  Notes: Light-time is already taken into account
 *  @param [in] Mjd_TDB         Modified julian date of TDB
 *  @param [out] r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun(geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF))
 */
void JPL_Eph_DE430(double Mjd_TBD, double *r_Mercury, double *r_Venus, double *r_Earth, double *r_Mars, double *r_Jupiter, double *r_Saturn, double *r_Uranus, double *r_Neptune, double *r_Pluto, double *r_Moon, double *r_Sun);

#endif/*JPL_EPH_DE430_H*/

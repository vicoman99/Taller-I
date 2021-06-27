/** @file gmst.h
 *  @brief Specification of gmst
 *
 *  @author Victor Coman
 *  @date May 2021
 */

#ifndef GMST_H
#define GMST_H

/** Greenwich Mean Sidereal Time
 *  @param [in] Mjd_UT1    Modified Julian Date UT1
 *  @return gmstime	   GMST in [rad]
 */
double gmst(double Mjd_UT1);

#endif/*GMST_H*/
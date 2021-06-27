/** @file timediff.h
 *  @brief Specification of timediff
 *
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef TIMEDIFF_H_
#define TIMEDIFF_H_

/** Time differences [s]
 *  @param [in] UT1_UTC, TAI_UTC
 *  @param [out] UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC
 */
void timediff(double UT1_UTC, double TAI_UTC, double *UT1_TAI, double *UTC_GPS, double *UT1_GPS, double *TT_UTC, double *GPS_UTC);

#endif
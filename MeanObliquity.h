/** @file MeanObliquity.h
 *  @brief Specification of MeanObliquity
 *
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef MEANOBLIQUITY_H
#define MEANOBLIQUITY_H

/** Computes the mean obliquity of the ecliptic
 *  @param [in] Mjd_TT    Modified Julian Date (Terrestrial Time)
 *  @return MOblq     Mean obliquity of the ecliptic [rad]
 */
double MeanObliquity(double Mjd_TT);

#endif

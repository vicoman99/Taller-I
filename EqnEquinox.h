/** @file EqnEquinox.h
 *  @brief Specification of EqnEquinox
 *
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef EQNEQUINOX_H
#define EQNEQUINOX_H

/** Computation of the equation of the equinoxes
 *  @brief Notes: The equation of the equinoxes dpsi*cos(eps) is the right ascension of the mean equinox referred to the true equator and equinox and is equal to the difference between apparent and mean sidereal time.
 *  @param [in] Mjd_TT    Modified Julian Date (Terrestrial Time)
 *  @return EqE      Equation of the equinoxes
 */
double EqnEquinox(double Mjd_TT);

#endif/*EQNEQUINOX_H*/
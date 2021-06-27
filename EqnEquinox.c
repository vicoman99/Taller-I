/** @file EqnEquinox.h
 *  @brief Implementation of EqnEquinox
 *
 *  @author Victor Coman
 *  @date May 2021
*/

#include "EqnEquinox.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include <math.h>

double EqnEquinox(double Mjd_TT)
{
    double dpsi, deps;
    // Nutation in longitude and obliquity
    NutAngles(Mjd_TT, &dpsi, &deps);

    // Equation of the equinoxes
    double EqE = dpsi * cos(MeanObliquity(Mjd_TT));

    return EqE;
}

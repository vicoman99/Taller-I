/** @file gast.h
 *  @brief Implementation of gast
 *  @author Victor Coman
 *  @date May 2021
*/
#include "gast.h"
#include "gmst.h"
#include <math.h>
#include "sat_const.h"
#include "EqnEquinox.h"

double gast(double Mjd_UT1)
{
    return fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2 * M_PI);
}

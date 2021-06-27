/** @file MeanObliquity.h
 *  @brief Implementation of MeanObliquity
 *
 *  @author Victor Coman
 *  @date May 2021
*/
#include <math.h>
#include "MeanObliquity.h"
#include "sat_const.h"

double MeanObliquity(double Mjd_TT)
{
    double T = (Mjd_TT - MJD_J2000) / 36525.0;
    double MOblq = Rad * (84381.448 / 3600.0 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0);

    return MOblq;
}

/** @file Mjday.h
 *  @brief Implementation of Mjday
 *  @author Victor Coman
 *  @date May 2021
*/
#include <math.h>
#include "Mjday.h"

double Mjday(int yr, int mon, int day, int hr, int min, int sec)
{
    double jd = 367.0 * yr - floor((7 * (yr + floor((mon + 9) / 12.0))) * 0.25) + floor(275 * mon / 9.0) + day + 1721013.5 + ((sec / 60.0 + min) / 60.0 + hr) / 24.0;
    double Mjd = jd - 2400000.5;
    return Mjd;
}

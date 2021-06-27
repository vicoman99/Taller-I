/** @file Mjday.h
 *  @brief Specification of Mjday
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef MJDAY_h_
#define MJDAY_h_

/** returns de julian date
 *  @param [in] year        - year
 *  @param [in] mon         - month
 *  @param [in] day         - day
 *  @param [in] hr          - universal time hour
 *  @param [in] min         - universal time min
 *  @param [in] sec         - universal time sec
 *  @return double
 */
double Mjday(int yr, int mon, int day, int hr, int min, int sec);

#endif

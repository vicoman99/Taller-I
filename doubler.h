/** @file doubler.h
 *  @brief Specification of doubler
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef DOUBLER_H
#define DOUBLER_H

/** Greenwich Apparent Sidereal Time
 *  @param [in] cc1
 *  @param [in] cc2
 *  @param [in] magrsite1
 *  @param [in] magrsite2
 *  @param [in] magr1in
 *  @param [in] magr2in
 *  @param [in] los1
 *  @param [in] los2
 *  @param [in] los3
 *  @param [in] rsite1
 *  @param [in] rsite2
 *  @param [in] rsite3
 *  @param [in] t1
 *  @param [in] t3
 *  @param [in] direct
 *  @param [in,out] r2
 *  @param [in,out] r3
 *  @param [in,out] f1
 *  @param [in,out] f2
 *  @param [in,out] magr1
 *  @param [in,out] magr2
 *  @param [in,out] a
 *  @param [in,out] deltae32
 */
void doubler(double cc1, double cc2, double magrsite1, double magrsite2, double magr1in, double magr2in, double *los1, double *los2, double *los3, double *rsite1, double *rsite2, double *rsite3, double t1, double t3, char direct, double *r2, double *r3, double *f1, double *f2, double *q1, double *magr1, double *magr2, double *a, double *deltae32);

#endif
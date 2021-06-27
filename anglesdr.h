/** @file anglesdr.h
 *  @brief Specification of anglesdr
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef ANGLESDR_H
#define ANGLESDR_H

/** this function solves the problem of orbit determination using three optical sightings.
 *  @param [in] az1      - azimuth at t1               rad
 *  @param [in] az2      - azimuth at t2               rad
 *  @param [in] az3      - azimuth at t3               rad
 *  @param [in] el1      - elevation at t1             rad
 *  @param [in] el2      - elevation at t2             rad
 *  @param [in] el3      - elevation at t3             rad
 *  @param [in] Mjd1     - Modified julian date of t1
 *  @param [in] Mjd2     - Modified julian date of t2
 *  @param [in] Mjd3     - Modified julian date of t3
 *  @param [in] rsite1   - ijk site1 position vector   m
 *  @param [in] rsite2   - ijk site2 position vector   m
 *  @param [in] rsite3   - ijk site3 position vector   m
 *  @param [out] r        - ijk position vector at t2   m
 *  @param [out] v        - ijk velocity vector at t2   m/s
 */
void anglesdr (double *r2, double *v2, double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, double *rsite1, double *rsite2, double *rsite3);

#endif

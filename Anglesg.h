/** @file anglesg.h
 *  @brief Specification of anglesg
 *
 *  @author Victor Coman
 *  @date May 2021
 */
#ifndef ANGLESG_H
#define ANLGESG_H

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
 *  @param [in] Rs1      - ijk site1 position vector   m
 *  @param [in] Rs2      - ijk site2 position vector   m
 *  @param [in] Rs3      - ijk site3 position vector   m
 *  @param [out] r        - ijk position vector at t2   m
 *  @param [out] v        - ijk velocity vector at t2   m/s
 */
void anglesg(double *r2, double *v2, double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, double *Rs1, double *Rs2, double *Rs3);

#endif /*ANGLESG_H*/
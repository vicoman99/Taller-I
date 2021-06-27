/** @file MeasUpdate.h
 *  @brief Implementation of MeasUpdate
 *  @author Victor Coman
 *  @date May 2021
*/

#include <stdio.h>
#include <math.h>
#include "arrays.h"
#include "sat_const.h"

void MeasUpdate(double *x, int nx, double z, double g, double s, double *G, int ng, double **P, int np, int n, double *K)
{
    double Inv_W = s; // Inverse weight (measurement covariance)

    // Kalman gain
    double parte1 = 1 / (Inv_W + dot(G, 6, matXvec(P, 6, 6, G, 6), 6));
    K = vec_x_esc(matXvec(P, 6, 6, G, 6), 6, parte1);
    //K = P*G'*inv(Inv_W+G*P*G');

    // State update
    x = sumV(x, 6, vec_x_esc(K, 6, z - g), 6);
    //x = x + K*(z-g);

    //Covariance update
    P = prod(sum(eye(6), 6, 6, mat_x_esc(vecXvecTrans(K, 6, G, 6), 6, -1), 6, 6), 6, 6, P, 6, 6);
    //P = (eye(n)-K*G)*P;
}

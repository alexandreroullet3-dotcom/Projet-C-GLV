#ifndef EC_ADD_PROJ_H
#define EC_ADD_PROJ_H

#include "EC_struct.h"

/*
 * =========================
 * Addition de deux points projectifs
 * 
 * Formules mathématiques (coordonnées Jacobiennes) :
 * 
 *   - Doublement : 
 *       X3 = M^2 - 2*S
 *       Y3 = M*(S - X3) - 8*Y^4
 *       Z3 = 2*Y*Z
 *     avec M = 3*X^2 + a*Z^4, S = 4*X*Y^2
 * 
 *   - Addition : 
 *       U1 = X1*Z2^2, U2 = X2*Z1^2
 *       S1 = Y1*Z2^3, S2 = Y2*Z1^3
 *       H = U2 - U1, R = S2 - S1
 *       X3 = R^2 - H^3 - 2*U1*H^2
 *       Y3 = R*(U1*H^2 - X3) - S1*H^3
 *       Z3 = H*Z1*Z2
 * =========================
 */

/* Doublement d'un point projectif : R = 2P */
void ec_point_double_proj(ECPointProj *R, const ECPointProj *P, const ECCurve *E);

/* Addition de deux points projectifs : R = P + Q */
void ec_point_add_proj(ECPointProj *R, const ECPointProj *P, const ECPointProj *Q, const ECCurve *E);

#endif /* EC_ADD_PROJ_H */

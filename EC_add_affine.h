#ifndef EC_ADD_AFFINE_H
#define EC_ADD_AFFINE_H

#include "EC_struct.h"

/*
 * =========================
 * Addition et doublement en coordonnées affines
 * =========================
 */

/* Addition de deux points affines : R = P + Q */
void ec_point_add_affine(ECPointAffine *R,
                         const ECPointAffine *P,
                         const ECPointAffine *Q,
                         const ECCurve *E);

/* Doublement d'un point affine : R = 2P */
void ec_point_double_affine(ECPointAffine *R, const ECPointAffine *P, const ECCurve *E);

/* Négation d'un point : R = -P */
void ec_point_affine_neg(ECPointAffine *R, const ECPointAffine *P, const ECCurve *E);

#endif /* EC_ADD_AFFINE_H */

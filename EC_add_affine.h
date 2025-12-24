#ifndef EC_ADD_AFFINE_H
#define EC_ADD_AFFINE_H 
#include "EC_struct.h"

void ec_point_add_affine(ECPointAffine *R, const ECPointAffine *P, const ECPointAffine *Q, const ECCurve *E);

void ec_point_double_affine(ECPointAffine *R, const ECPointAffine *P, const ECCurve *E);
#endif
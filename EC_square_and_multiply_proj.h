#ifndef EC_SQUARE_AND_MULTIPLY_PROJ_H
#define EC_SQUARE_AND_MULTIPLY_PROJ_H

#include "EC_add_proj.h"

// multiplication scalaire : R = k * P
// k est un mpz_t
void ec_scalar_mul_proj(ECPointProj *R, const ECPointProj *P, const mpz_t k, const ECCurve *E);

#endif
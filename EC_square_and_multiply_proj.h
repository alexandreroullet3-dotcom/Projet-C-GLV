#ifndef EC_SQUARE_AND_MULTIPLY_PROJ_H
#define EC_SQUARE_AND_MULTIPLY_PROJ_H

#include "EC_add_proj.h"

/*
 * =========================
 * Multiplication scalaire sur points projectifs
 * Algorithme square-and-multiply classique
 * =========================
 *
 * R = k * P, où k est un mpz_t
 */
void ec_scalar_mul_proj(ECPointProj *R,
                        const ECPointProj *P,
                        const mpz_t k,
                        const ECCurve *E);

#endif /* EC_SQUARE_AND_MULTIPLY_PROJ_H */

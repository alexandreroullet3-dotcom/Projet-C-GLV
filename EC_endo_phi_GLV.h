#ifndef EC_ENDO_PHI_GLV_H
#define EC_ENDO_PHI_GLV_H

#include "EC_struct.h"

void ec_endo_phi_point_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t beta);

void ec_glv_decompose(mpz_t k1, mpz_t k2, const mpz_t k, const mpz_t n, const mpz_t lambda,
                    const mpz_t x1, const mpz_t y1, const mpz_t x2, const mpz_t y2);

#endif
#ifndef EC_GLV_H
#define EC_GLV_H

#include "EC_endo_phi_GLV.h"
#include "double_scalar_multiplication.h"
#include "EC_conversions.h"
#include "glv_decompose.h"

void ec_scal_mul_glv(ECPointProj *R, const ECPointProj *P, const mpz_t k, const ECCurve *E, 
    const Z2 *v1, const Z2 *v2, const mpz_t beta, mpz_t n, int type);
#endif
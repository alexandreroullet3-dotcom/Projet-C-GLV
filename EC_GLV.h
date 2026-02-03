#ifndef EC_GLV_H
#define EC_GLV_H

#include "EC_endo_phi_GLV.h"
#include "double_scalar_multiplication.h"
#include "glv_decompose.h"

void ec_scal_mul_glv(ECPointProj *R, const ECPointProj *P, const ECPointProj *phiP, const mpz_t k, const ECCurve *E, 
    const Z2 *v1, const Z2 *v2);
#endif
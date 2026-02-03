#include "EC_GLV.h"

// Multiplication scalaire GLV  R = k * P

void ec_scal_mul_glv(ECPointProj *R, const ECPointProj *P, const ECPointProj *phiP, const mpz_t k, const ECCurve *E, 
    const Z2 *v1, const Z2 *v2){
    // on commence par decomposer k
    Z2 v;
    z2_init(&v);
    glv_nearest_vector(&v, k, v1, v2);

    ec_double_scalar_multiplication(R, P, phiP, v.x, v.y, 2, E);
    
    z2_clear(&v);
}
#include "EC_GLV.h"

// Multiplication scalaire GLV  R = k * P

void ec_scal_mul_glv(ECPointProj *R, const ECPointProj *P, const mpz_t k, const ECCurve *E, 
    const Z2 *v1, const Z2 *v2, const mpz_t beta, mpz_t n){
    // on commence par decomposer k
    Z2 v;
    z2_init(&v);
    glv_nearest_vector(&v, k, v1, v2, n);

    // Calculons Q = phi(P) 
    //pour cela on passe en affine pour simplifier les calcules et appliquer seulement x*beta
    ECPointAffine P_aff, Q_aff;
    ECPointProj Q;
    ec_point_affine_init(&P_aff);
    ec_point_affine_init(&Q_aff);
    ec_point_proj_init(&Q);
    proj_to_affine(&P_aff, P, E);
    ec_endo_phi_point_affine(&Q_aff, &P_aff, E, beta);
    affine_to_proj(&Q, &Q_aff);


    ec_double_scalar_multiplication(R, P, &Q, v.x, v.y, 1, E);
    
    z2_clear(&v);
    ec_point_proj_clear(&Q);
    ec_point_affine_clear(&P_aff);
    ec_point_affine_clear(&Q_aff);
}
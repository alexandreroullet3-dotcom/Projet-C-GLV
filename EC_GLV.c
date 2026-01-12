#include <stdlib.h>
#include <gmp.h>
#include "EC_struct.h"
#include "EC_add_proj.h"
#include "EC_endo_phi_GLV.h"
#include "double_scalar_multiplication.h"
#include "EC_conversions.h"

// Multiplication scalaire GLV  R = k * P

void ec_scal_mul_glv(ECPointProj *R, const ECPointProj *P, const mpz_t k, const ECCurve * E, const mpz_t x1, 
                    const mpz_t y1, const mpz_t x2, const mpz_t y2, const mpz_t beta){
    // on commence par decomposer k
    mpz_t k1,k2;
    mpz_inits(k1,k2,NULL);
    ec_glv_decompose(k1,k2,k,x1,y1,x2,y2);

    // Calculons Q = phi(P) 
    //pour cela on passe en affine pour simplifier les calcules et appliquer seulement x*beta
    ECPointAffine P_aff, Q_aff;
    ECPointProj Q;
    ec_point_affine_init(&P_aff);
    ec_point_affine_init(&Q_aff);
    ec_point_proj_init(&Q);
    proj_to_affine(&P_aff, P,E);
    ec_endo_phi_point_affine(&Q_aff,&P_aff,E,beta);
    affine_to_proj(&Q, &Q_aff);


    ec_double_scalar_multiplication(R, P, &Q, k1, k2, 2, E);

    mpz_clears(k1, k2, NULL );
    ec_point_proj_clear(&Q);
    ec_point_affine_clear(&P_aff);
    ec_point_affine_clear(&Q_aff);
}

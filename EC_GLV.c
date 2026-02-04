#include "EC_GLV.h"

/*
 * Module GLV: multiplication scalaire avec endomorphisme.
 */
void ec_scal_mul_glv(ECPointProj *R,
                     const ECPointProj *P,
                     const ECPointProj *phiP,
                     const mpz_t k,
                     const ECCurve *E,
                     const Z2 *v1,
                     const Z2 *v2)
{
    // Décomposition de k dans la base courte (v1, v2).
    Z2 v;
    z2_init(&v);
    glv_nearest_vector(&v, k, v1, v2);

    // Double multiplication scalaire via l'endomorphisme.
    // On a choisi ici w=2, cest ce qui donne empiriquement la meilleure accélération
    ec_double_scalar_multiplication(R, P, phiP, v.x, v.y, 2, E);

    z2_clear(&v);
}
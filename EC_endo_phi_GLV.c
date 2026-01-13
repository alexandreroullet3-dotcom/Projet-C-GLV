#include "EC_endo_phi_GLV.h"


// fonction qui applique l'endomorphisme phi sur un point P
// il est défini par phi(x,y) = (beta*x mod p, y)  avec beta racine cubiuqe de l'unité
void ec_endo_phi_point_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t beta){
    if(P->infinity){
        Q->infinity =1;
    }else{
    // Q.x =beta *P.x mod p
    mpz_mul(Q->x, P->x, beta);
    mpz_mod(Q->x, Q->x, E->p);

    // Q.y = P.y
    mpz_set(Q->y, P->y);
    Q-> infinity =0;
    }
}
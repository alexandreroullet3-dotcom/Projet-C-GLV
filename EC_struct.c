#include <gmp.h>
#include "EC_struct.h"


void ec_point_affine_init(ECPointAffine *P) {
    mpz_init(P->x);
    mpz_init(P->y);
    P->infinity = 1;
}

void ec_point_affine_clear(ECPointAffine *P) {
    mpz_clear(P->x);
    mpz_clear(P->y);
}

void ec_point_proj_init(ECPointProj *P) {
    mpz_init(P->X);
    mpz_init(P->Y);
    mpz_init(P->Z);
    P->infinity = 1;
}

void ec_point_proj_clear(ECPointProj *P) {
    mpz_clear(P->X);
    mpz_clear(P->Y);
    mpz_clear(P->Z);
}

void ec_curve_init(ECCurve *E) {
    mpz_init(E->p);
    mpz_init(E->a);
    mpz_init(E->b);
}

void ec_curve_clear(ECCurve *E) {
    mpz_clear(E->p);
    mpz_clear(E->a);
    mpz_clear(E->b);
}
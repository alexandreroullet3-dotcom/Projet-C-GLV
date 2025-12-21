#include <gmp.h>
#include "EC_struct.h"


void ec_point_init(ECPoint *P) {
    mpz_init(P->x);
    mpz_init(P->y);
    P->infinity = 1;
}

void ec_point_clear(ECPoint *P) {
    mpz_clear(P->x);
    mpz_clear(P->y);
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
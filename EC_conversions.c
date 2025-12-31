#include "EC_conversions.h"

void affine_to_proj(ECPointProj *R, const ECPointAffine *P)
{
    if (P->infinity) {
        R->infinity = 1;
        return;
    }
    mpz_set(R->X, P->x);
    mpz_set(R->Y, P->y);
    mpz_set_ui(R->Z, 1);
    R->infinity = 0;
}

void proj_to_affine(ECPointAffine *R,
                    const ECPointProj *P,
                    const ECCurve *E)
{
    if (P->infinity) {
        R->infinity = 1;
        return;
    }

    mpz_t Zinv, Z2, Z3;
    mpz_inits(Zinv, Z2, Z3, NULL);

    mpz_invert(Zinv, P->Z, E->p);
    mpz_mul(Z2, Zinv, Zinv);
    mpz_mod(Z2, Z2, E->p);

    mpz_mul(Z3, Z2, Zinv);
    mpz_mod(Z3, Z3, E->p);

    mpz_mul(R->x, P->X, Z2);
    mpz_mod(R->x, R->x, E->p);

    mpz_mul(R->y, P->Y, Z3);
    mpz_mod(R->y, R->y, E->p);

    R->infinity = 0;

    mpz_clears(Zinv, Z2, Z3, NULL);
}
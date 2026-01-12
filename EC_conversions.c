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

int ec_cmp_aff(ECPointAffine *P, ECPointAffine *Q){
    if (!(mpz_cmp(P->x, Q->x)||mpz_cmp(P->y,Q->y)))
        return 1;
    else
        return 0;
}

int ec_cmp_proj(ECPointProj *P, ECPointProj *Q, ECCurve *E){
    ECPointAffine Paff, Qaff;
    ec_point_affine_init(&Paff);
    ec_point_affine_init(&Qaff);
    proj_to_affine(&Paff, P, E);
    proj_to_affine(&Qaff, Q, E);
    int i = ec_cmp_aff(&Paff, &Qaff);
    ec_point_affine_clear(&Paff);
    ec_point_affine_clear(&Qaff);
    return i; 
}
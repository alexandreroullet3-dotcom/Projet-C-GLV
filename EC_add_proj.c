#include "EC_add_proj.h"

/*
 * =========================
 * Doublement d'un point projectif
 * Formules Jacobiennes pour courbes y^2 = x^3 + a2*x^2 + a4*x + a6
 * =========================
 */
void ec_point_double_proj(ECPointProj *R,
                          const ECPointProj *P,
                          const ECCurve *E)
{
    if (P->infinity || mpz_cmp_ui(P->Y, 0) == 0) {
        R->infinity = 1;
        mpz_set_ui(R->Z, 0);
        return;
    }

    mpz_t Y2, S, M, Z2, Z4, tmp;
    mpz_inits(Y2, S, M, Z2, Z4, tmp, NULL);

    /* Y2 = Y^2 */
    mpz_mul(Y2, P->Y, P->Y);
    mpz_mod(Y2, Y2, E->p);

    /* S = 4 X Y^2 */
    mpz_mul(S, P->X, Y2);
    mpz_mul_ui(S, S, 4);
    mpz_mod(S, S, E->p);

    /* Z2 = Z^2 */
    mpz_mul(Z2, P->Z, P->Z);

    /* Z4 = Z^4 */
    mpz_mul(Z4, Z2, Z2);

    /* M = 3X^2 */
    mpz_mul(M, P->X, P->X);
    mpz_mul_ui(M, M, 3);

    /* + 2 a2 X Z^2 */
    mpz_mul(tmp, E->a2, P->X);
    mpz_mul(tmp, tmp, Z2);
    mpz_mul_ui(tmp, tmp, 2);
    mpz_add(M, M, tmp);

    /* + a Z^4 */
    mpz_mul(tmp, E->a, Z4);
    mpz_add(M, M, tmp);
    mpz_mod(M, M, E->p);

    /* Z3 = 2 Y Z */
    mpz_mul(R->Z, P->Y, P->Z);
    mpz_mul_ui(R->Z, R->Z, 2);
    mpz_mod(R->Z, R->Z, E->p);

    /* X3 = M^2 − 2S − a2 Z3^2 */
    mpz_mul(R->X, M, M);

    mpz_sub(R->X, R->X, S);
    mpz_sub(R->X, R->X, S);

    mpz_mul(tmp, R->Z, R->Z);
    mpz_mul(tmp, tmp, E->a2);
    mpz_sub(R->X, R->X, tmp);
    mpz_mod(R->X, R->X, E->p);

    /* Y3 = M(S − X3) − 8Y^4 */
    mpz_sub(tmp, S, R->X);
    mpz_mul(tmp, tmp, M);

    mpz_mul(Y2, Y2, Y2);      // Y^4
    mpz_mul_ui(Y2, Y2, 8);

    mpz_sub(R->Y, tmp, Y2);
    mpz_mod(R->Y, R->Y, E->p);

    R->infinity = 0;

    mpz_clears(Y2, S, M, Z2, Z4, tmp, NULL);
}


/*
 * =========================
 * Addition de deux points projectifs
 * Formules Jacobiennes pour courbes y^2 = x^3 + a2*x^2 + a4*x + a6
 * =========================
 */
void ec_point_add_proj(ECPointProj *R,
                       const ECPointProj *P,
                       const ECPointProj *Q,
                       const ECCurve *E)
{
    if (P->infinity) {
        ec_point_proj_copy(R, Q);
        return;
    }
    if (Q->infinity) {
        ec_point_proj_copy(R, P);
        return;
    }

    mpz_t Z1Z1, Z2Z2, U1, U2, S1, S2;
    mpz_t H, Rr, H2, H3, U1H2, tmp;

    mpz_inits(Z1Z1, Z2Z2, U1, U2, S1, S2,
              H, Rr, H2, H3, U1H2, tmp, NULL);

    /* Z1^2, Z2^2 */
    mpz_mul(Z1Z1, P->Z, P->Z);
    mpz_mul(Z2Z2, Q->Z, Q->Z);

    /* U1, U2 */
    mpz_mul(U1, P->X, Z2Z2);
    mpz_mul(U2, Q->X, Z1Z1);

    /* S1 = Y1 Z2^3 */
    mpz_mul(S1, Z2Z2, Q->Z);
    mpz_mul(S1, S1, P->Y);

    /* S2 = Y2 Z1^3 */
    mpz_mul(S2, Z1Z1, P->Z);
    mpz_mul(S2, S2, Q->Y);

    mpz_mod(U1, U1, E->p);
    mpz_mod(U2, U2, E->p);
    mpz_mod(S1, S1, E->p);
    mpz_mod(S2, S2, E->p);

    /* cas U1 = U2 */
    if (mpz_cmp(U1, U2) == 0) {
        if (mpz_cmp(S1, S2) != 0) {
            R->infinity = 1;
            mpz_set_ui(R->Z, 0);
        } else {
            ec_point_double_proj(R, P, E);
        }
        goto clear;
    }

    /* H = U2 − U1 */
    mpz_sub(H, U2, U1);

    /* R = S2 − S1 */
    mpz_sub(Rr, S2, S1);

    mpz_mod(H, H, E->p);
    mpz_mod(Rr, Rr, E->p);

    /* H², H³ */
    mpz_mul(H2, H, H);
    mpz_mul(H3, H2, H);

    /* U1H² */
    mpz_mul(U1H2, U1, H2);

    /* Z3 = H Z1 Z2 */
    mpz_mul(R->Z, H, P->Z);
    mpz_mul(R->Z, R->Z, Q->Z);
    mpz_mod(R->Z, R->Z, E->p);

    /* X3 */
    mpz_mul(R->X, Rr, Rr);
    mpz_sub(R->X, R->X, H3);
    mpz_sub(R->X, R->X, U1H2);
    mpz_sub(R->X, R->X, U1H2);

    mpz_mul(tmp, R->Z, R->Z);
    mpz_mul(tmp, tmp, E->a2);
    mpz_sub(R->X, R->X, tmp);
    mpz_mod(R->X, R->X, E->p);

    /* Y3 */
    mpz_sub(tmp, U1H2, R->X);
    mpz_mul(tmp, tmp, Rr);

    mpz_mul(S1, S1, H3);
    mpz_sub(R->Y, tmp, S1);
    mpz_mod(R->Y, R->Y, E->p);

    R->infinity = 0;

clear:
    mpz_clears(Z1Z1, Z2Z2, U1, U2, S1, S2,
               H, Rr, H2, H3, U1H2, tmp, NULL);
}


/*
 * Négation d'un point projectif.
 */
void ec_point_proj_neg(ECPointProj *R, const ECPointProj *P)
{
    if (P->infinity) {
        R->infinity = 1;
    } else {
        R->infinity = 0;
        mpz_set(R->X, P->X);
        mpz_neg(R->Y, P->Y);
        mpz_set(R->Z, P->Z);
    }
}
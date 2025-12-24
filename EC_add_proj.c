#include "EC_add_proj.h"

void ec_point_double_proj(ECPointProj *R,
                     const ECPointProj *P,
                     const ECCurve *E)
{
    if (P->infinity || mpz_cmp_ui(P->Y, 0) == 0) {
        R->infinity = 1;
        return;
    }

    mpz_t S, M, T;
    mpz_inits(S, M, T, NULL);

    // S = 4 * X * Y^2
    mpz_mul(S, P->Y, P->Y);
    mpz_mul(S, S, P->X);
    mpz_mul_ui(S, S, 4);
    mpz_mod(S, S, E->p);

    // M = 3*X^2 + a*Z^4
    mpz_mul(M, P->X, P->X);
    mpz_mul_ui(M, M, 3);

    mpz_mul(T, P->Z, P->Z);
    mpz_mul(T, T, T);
    mpz_mul(T, T, E->a);

    mpz_add(M, M, T);
    mpz_mod(M, M, E->p);

    // X3 = M^2 - 2*S
    mpz_mul(R->X, M, M);
    mpz_sub(R->X, R->X, S);
    mpz_sub(R->X, R->X, S);
    mpz_mod(R->X, R->X, E->p);

    // Y3 = M*(S - X3) - 8*Y^4
    mpz_sub(R->Y, S, R->X);
    mpz_mul(R->Y, M, R->Y);

    mpz_mul(T, P->Y, P->Y);
    mpz_mul(T, T, T);
    mpz_mul_ui(T, T, 8);

    mpz_sub(R->Y, R->Y, T);
    mpz_mod(R->Y, R->Y, E->p);

    // Z3 = 2*Y*Z
    mpz_mul(R->Z, P->Y, P->Z);
    mpz_mul_ui(R->Z, R->Z, 2);
    mpz_mod(R->Z, R->Z, E->p);

    R->infinity = 0;

    mpz_clears(S, M, T, NULL);
}

void ec_point_add_proj(ECPointProj *R,
                  const ECPointProj *P,
                  const ECPointProj *Q,
                  const ECCurve *E)
{
    if (P->infinity) { *R = *Q; return; }
    if (Q->infinity) { *R = *P; return; }

    mpz_t U1, U2, S1, S2, H, Rr;
    mpz_inits(U1, U2, S1, S2, H, Rr, NULL);

    // U1 = X1*Z2^2
    mpz_mul(U1, Q->Z, Q->Z);
    mpz_mul(U1, U1, P->X);
    mpz_mod(U1, U1, E->p);

    // U2 = X2*Z1^2
    mpz_mul(U2, P->Z, P->Z);
    mpz_mul(U2, U2, Q->X);
    mpz_mod(U2, U2, E->p);

    // S1 = Y1*Z2^3
    mpz_mul(S1, Q->Z, Q->Z);
    mpz_mul(S1, S1, Q->Z);
    mpz_mul(S1, S1, P->Y);
    mpz_mod(S1, S1, E->p);

    // S2 = Y2*Z1^3
    mpz_mul(S2, P->Z, P->Z);
    mpz_mul(S2, S2, P->Z);
    mpz_mul(S2, S2, Q->Y);
    mpz_mod(S2, S2, E->p);

    if (mpz_cmp(U1, U2) == 0) {
        if (mpz_cmp(S1, S2) != 0) {
            R->infinity = 1;
        } else {
            ec_point_double_proj(R, P, E);
        }
        mpz_clears(U1, U2, S1, S2, H, Rr, NULL);
        return;
    }

    mpz_sub(H, U2, U1);
    mpz_sub(Rr, S2, S1);

    // X3 = Rr^2 - H^3 - 2*U1*H^2
    mpz_mul(R->X, Rr, Rr);
    mpz_mul(H, H, H);
    mpz_sub(R->X, R->X, H);
    mpz_sub(R->X, R->X, H);
    mpz_mod(R->X, R->X, E->p);

    // Y3 = Rr*(U1*H^2 - X3) - S1*H^3
    mpz_sub(R->Y, H, R->X);
    mpz_mul(R->Y, Rr, R->Y);
    mpz_mod(R->Y, R->Y, E->p);

    // Z3 = H*Z1*Z2
    mpz_mul(R->Z, H, P->Z);
    mpz_mul(R->Z, R->Z, Q->Z);
    mpz_mod(R->Z, R->Z, E->p);

    R->infinity = 0;

    mpz_clears(U1, U2, S1, S2, H, Rr, NULL);
}





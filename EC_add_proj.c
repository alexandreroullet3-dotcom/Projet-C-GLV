#include "EC_add_proj.h"

/*
 * =========================
 * Doublement d'un point projectif
 * Formules Jacobiennes pour courbes y^2 = x^3 + a2*x^2 + a4*x + a6
 * =========================
 */
void ec_point_double_proj(ECPointProj *R, const ECPointProj *P, const ECCurve *E)
{
    if (P->infinity || mpz_cmp_ui(P->Y, 0) == 0) {
        R->infinity = 1;
        mpz_set_ui(R->Z, 0);
        return;
    }

    mpz_t T1, T2, M, S, Ta2;
    mpz_t X3, Y3, Z3;
    mpz_inits(T1, T2, M, S, Ta2, X3, Y3, Z3, NULL);

    // 1. Calcul de T1 = Y^2 et S = 4*X*Y^2
    mpz_mul(T1, P->Y, P->Y);
    mpz_mod(T1, T1, E->p);
    mpz_mul(S, P->X, T1);
    mpz_mul_ui(S, S, 4);
    mpz_mod(S, S, E->p);

    // 2. Calcul de M = 3*X^2 + 2*a2*X*Z^2 + a*Z^4
    mpz_mul(M, P->X, P->X);
    mpz_mul_ui(M, M, 3);      // 3*X^2
    mpz_mul(T2, P->Z, P->Z);  // Z^2
    
    // Terme 2*a2*X*Z^2
    mpz_mul(Ta2, E->a2, P->X);
    mpz_mul(Ta2, Ta2, T2);
    mpz_mul_ui(Ta2, Ta2, 2);
    mpz_add(M, M, Ta2);       // + 2*a2*X*Z^2
    
    // Terme a*Z^4
    mpz_mul(T2, T2, T2);      // Z^4
    mpz_mul(T2, T2, E->a);    // a*Z^4
    mpz_add(M, M, T2);        // + a*Z^4
    mpz_mod(M, M, E->p);

    // 3. Calcul de Z3 = 2*Y*Z (Nécessaire pour corriger X3)
    mpz_mul(Z3, P->Y, P->Z);
    mpz_mul_ui(Z3, Z3, 2);
    mpz_mod(Z3, Z3, E->p);

    // 4. Calcul de X3 = M^2 - 2*S - a2*Z3^2
    mpz_mul(X3, M, M);
    mpz_sub(X3, X3, S);
    mpz_sub(X3, X3, S);
    
    // Correction a2
    mpz_mul(Ta2, Z3, Z3);
    mpz_mul(Ta2, Ta2, E->a2);
    mpz_sub(X3, X3, Ta2);
    mpz_mod(X3, X3, E->p);

    // 5. Calcul de Y3 = M*(S - X3) - 8*Y^4
    mpz_mul(T2, T1, T1);      // Y^4
    mpz_mul_ui(T2, T2, 8);    // 8*Y^4
    mpz_sub(Y3, S, X3);
    mpz_mul(Y3, M, Y3);
    mpz_sub(Y3, Y3, T2);
    mpz_mod(Y3, Y3, E->p);

    mpz_set(R->X, X3);
    mpz_set(R->Y, Y3);
    mpz_set(R->Z, Z3);
    R->infinity = 0;

    mpz_clears(T1, T2, M, S, Ta2, X3, Y3, Z3, NULL);
}

/*
 * =========================
 * Addition de deux points projectifs
 * Formules Jacobiennes pour courbes y^2 = x^3 + a2*x^2 + a4*x + a6
 * =========================
 */
void ec_point_add_proj(ECPointProj *R, const ECPointProj *P, 
                       const ECPointProj *Q, const ECCurve *E)
{
    if (P->infinity) { ec_point_proj_copy(R, Q); return; }
    if (Q->infinity) { ec_point_proj_copy(R, P); return; }

    mpz_t U1, U2, S1, S2, H, Rr, H2, H3, U1H2, Ta2;
    mpz_t X3, Y3, Z3;
    mpz_inits(U1, U2, S1, S2, H, Rr, H2, H3, U1H2, Ta2, X3, Y3, Z3, NULL);

    // 1. Calcul des coordonnées intermédiaires
    mpz_mul(H2, Q->Z, Q->Z); mpz_mod(H2, H2, E->p);
    mpz_mul(U1, P->X, H2);   mpz_mod(U1, U1, E->p);
    mpz_mul(S1, P->Y, H2);   mpz_mul(S1, S1, Q->Z); mpz_mod(S1, S1, E->p);

    mpz_mul(H2, P->Z, P->Z); mpz_mod(H2, H2, E->p);
    mpz_mul(U2, Q->X, H2);   mpz_mod(U2, U2, E->p);
    mpz_mul(S2, Q->Y, H2);   mpz_mul(S2, S2, P->Z); mpz_mod(S2, S2, E->p);

    if (mpz_cmp(U1, U2) == 0) {
        if (mpz_cmp(S1, S2) != 0) {
            R->infinity = 1;
            mpz_set_ui(R->Z, 0);
        } else {
            ec_point_double_proj(R, P, E);
        }
        mpz_clears(U1, U2, S1, S2, H, Rr, H2, H3, U1H2, Ta2, X3, Y3, Z3, NULL);
        return;
    }

    mpz_sub(H, U2, U1);
    mpz_sub(Rr, S2, S1);
    mpz_mul(H2, H, H);       mpz_mod(H2, H2, E->p);
    mpz_mul(H3, H2, H);      mpz_mod(H3, H3, E->p);
    mpz_mul(U1H2, U1, H2);   mpz_mod(U1H2, U1H2, E->p);

    // 2. Calcul de Z3 D'ABORD (Nécessaire pour Ta2)
    mpz_mul(Z3, H, P->Z);
    mpz_mul(Z3, Z3, Q->Z);
    mpz_mod(Z3, Z3, E->p);

    // 3. Calcul de X3 avec le terme correctif a2
    mpz_mul(X3, Rr, Rr);
    mpz_sub(X3, X3, H3);
    mpz_sub(X3, X3, U1H2);
    mpz_sub(X3, X3, U1H2);
    
    // Correction pour courbe généralisée : - a2 * Z3^2
    mpz_mul(Ta2, Z3, Z3);
    mpz_mul(Ta2, Ta2, E->a2);
    mpz_sub(X3, X3, Ta2);
    mpz_mod(X3, X3, E->p);

    // 4. Calcul de Y3 utilisant le X3 DÉJÀ CORRIGÉ
    mpz_sub(Y3, U1H2, X3);
    mpz_mul(Y3, Rr, Y3);
    mpz_mul(S1, S1, H3);
    mpz_sub(Y3, Y3, S1);
    mpz_mod(Y3, Y3, E->p);

    mpz_set(R->X, X3);
    mpz_set(R->Y, Y3);
    mpz_set(R->Z, Z3);
    R->infinity = 0;

    mpz_clears(U1, U2, S1, S2, H, Rr, H2, H3, U1H2, Ta2, X3, Y3, Z3, NULL);
}

void ec_point_proj_neg(ECPointProj *R, const ECPointProj *P){
    if (P->infinity){
        R->infinity = 1;
    }
    else{
        R->infinity = 0;
        mpz_set(R->X, P->X);
        mpz_neg(R->Y, P->Y);
        mpz_set(R->Z, P->Z);
    }
}
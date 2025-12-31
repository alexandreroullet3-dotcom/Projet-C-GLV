#include "EC_add_proj.h"


void ec_point_double_proj(ECPointProj *R, const ECPointProj *P, const ECCurve *E)
{
    // 1. Gestion du point à l'infini ou Y=0
    if (P->infinity || mpz_cmp_ui(P->Y, 0) == 0) {
        R->infinity = 1;
        mpz_set_ui(R->Z, 0);
        return;
    }

    // 2. Initialisation des variables temporaires
    // On a besoin de X3, Y3, Z3 pour ne pas écrire dans R trop tôt et donc le remplacé 
    mpz_t T1, T2, M, S; 
    mpz_t X3, Y3, Z3;
    mpz_inits(T1, T2, M, S, X3, Y3, Z3, NULL);

    // Formules standard (Jacobiennes) :
    // A = X^2
    // B = Y^2
    // C = B^2
    // S = 4 * X * B
    // M = 3 * A + a * Z^4
    
    // T1 = Y^2
    mpz_mul(T1, P->Y, P->Y);
    mpz_mod(T1, T1, E->p);

    // S = 4 * X * Y^2
    mpz_mul(S, P->X, T1);
    mpz_mul_ui(S, S, 4);
    mpz_mod(S, S, E->p);

    // M = 3 * X^2 + a * Z^4
    mpz_mul(M, P->X, P->X);   // X^2
    mpz_mul_ui(M, M, 3);      // 3*X^2
    
    mpz_mul(T2, P->Z, P->Z);  // Z^2
    mpz_mul(T2, T2, T2);      // Z^4
    mpz_mul(T2, T2, E->a);    // a*Z^4
    
    mpz_add(M, M, T2);        
    mpz_mod(M, M, E->p);      // M total

    // Calcul des Résultats dans les variables temporaires 

    // X3 = M^2 - 2*S
    mpz_mul(X3, M, M);
    mpz_sub(X3, X3, S);
    mpz_sub(X3, X3, S);
    mpz_mod(X3, X3, E->p);

    // Z3 = 2 * Y * Z
    mpz_mul(Z3, P->Y, P->Z);
    mpz_mul_ui(Z3, Z3, 2);
    mpz_mod(Z3, Z3, E->p);

    // Y3 = M(S - X3) - 8*Y^4
    // T2 = Y^4 (T1 au carré)
    mpz_mul(T2, T1, T1);
    mpz_mul_ui(T2, T2, 8); // 8*Y^4

    mpz_sub(Y3, S, X3);
    mpz_mul(Y3, M, Y3);
    mpz_sub(Y3, Y3, T2);
    mpz_mod(Y3, Y3, E->p);

    // 3. Copie finale (C'est seulement ICI qu'on touche à R)
    mpz_set(R->X, X3);
    mpz_set(R->Y, Y3);
    mpz_set(R->Z, Z3);
    R->infinity = 0;

    // 4. Nettoyage
    mpz_clears(T1, T2, M, S, X3, Y3, Z3, NULL);
}


void ec_point_add_proj(ECPointProj *R, const ECPointProj *P, const ECPointProj *Q, const ECCurve *E)
{
    if (P->infinity) { ec_point_proj_copy(R, Q); return; }
    if (Q->infinity) { ec_point_proj_copy(R, P); return; }

    mpz_t U1, U2, S1, S2, H, Rr;
    mpz_t X3, Y3, Z3; // Variables de résultat temporaires
    mpz_inits(U1, U2, S1, S2, H, Rr, X3, Y3, Z3, NULL);

    // U1 = X1 * Z2^2
    mpz_mul(U1, Q->Z, Q->Z);
    mpz_mod(U1, U1, E->p);
    mpz_mul(U1, U1, P->X);
    mpz_mod(U1, U1, E->p);

    // U2 = X2 * Z1^2
    mpz_mul(U2, P->Z, P->Z);
    mpz_mod(U2, U2, E->p);
    mpz_mul(U2, U2, Q->X);
    mpz_mod(U2, U2, E->p);

    // S1 = Y1 * Z2^3
    mpz_mul(S1, Q->Z, Q->Z);
    mpz_mul(S1, S1, Q->Z);
    mpz_mul(S1, S1, P->Y);
    mpz_mod(S1, S1, E->p);

    // S2 = Y2 * Z1^3
    mpz_mul(S2, P->Z, P->Z);
    mpz_mul(S2, S2, P->Z);
    mpz_mul(S2, S2, Q->Y);
    mpz_mod(S2, S2, E->p);

    if (mpz_cmp(U1, U2) == 0) {
        if (mpz_cmp(S1, S2) != 0) {
            R->infinity = 1;
            mpz_set_ui(R->Z, 0);
        } else {
            ec_point_double_proj(R, P, E);
        }
        mpz_clears(U1, U2, S1, S2, H, Rr, X3, Y3, Z3, NULL);
        return;
    }

    mpz_sub(H, U2, U1);      // H = U2 - U1
    mpz_sub(Rr, S2, S1);     // Rr = S2 - S1

    // Calcul avec variables temporaires pour H^2, H^3
    mpz_t H2, H3, U1H2;
    mpz_inits(H2, H3, U1H2, NULL);

    mpz_mul(H2, H, H);       
    mpz_mod(H2, H2, E->p);   // H^2

    mpz_mul(H3, H2, H);      
    mpz_mod(H3, H3, E->p);   // H^3

    mpz_mul(U1H2, U1, H2);
    mpz_mod(U1H2, U1H2, E->p); // U1*H^2

    // Calcul de X3 (dans temp)
    mpz_mul(X3, Rr, Rr);
    mpz_sub(X3, X3, H3);
    mpz_sub(X3, X3, U1H2);
    mpz_sub(X3, X3, U1H2);
    mpz_mod(X3, X3, E->p);

    // Calcul de Y3 (dans temp)
    mpz_sub(Y3, U1H2, X3);
    mpz_mul(Y3, Rr, Y3);
    mpz_mul(S1, S1, H3);     // S1*H^3
    mpz_sub(Y3, Y3, S1);
    mpz_mod(Y3, Y3, E->p);

    // Calcul de Z3 (dans temp)
    mpz_mul(Z3, H, P->Z);
    mpz_mul(Z3, Z3, Q->Z);
    mpz_mod(Z3, Z3, E->p);

    // --- Copie Finale ---
    mpz_set(R->X, X3);
    mpz_set(R->Y, Y3);
    mpz_set(R->Z, Z3);
    R->infinity = 0;

    mpz_clears(U1, U2, S1, S2, H, Rr, NULL);
    mpz_clears(H2, H3, U1H2, NULL);
    mpz_clears(X3, Y3, Z3, NULL);
}

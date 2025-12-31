#include "precompute_table.h"
#include "EC_add_proj.h"

/**
 * Pré-calcul de la table T[i][j] = i*P + j*Q
 * Retourne le tableau T[i][j], et éventuellement Ptab/Qtab si besoin
 */
ECPointProj** precompute_table(const ECPointProj *P,
                                      const ECPointProj *Q,
                                      unsigned int w,
                                      const ECCurve *E)
{
    unsigned int M = 1U << w;

    // --- Allocation et init de Ptab et Qtab ---
    ECPointProj *Ptab = malloc(M * sizeof(ECPointProj));
    ECPointProj *Qtab = malloc(M * sizeof(ECPointProj));
    for (unsigned int i = 0; i < M; i++) {
        ec_point_proj_init(&Ptab[i]);
        ec_point_proj_init(&Qtab[i]);
    }
    mpz_set_ui(Ptab[0].X, 0);
    mpz_set_ui(Ptab[0].Y, 1);
    mpz_set_ui(Ptab[0].Z, 0);
    Ptab[0].infinity = 1;
    mpz_set_ui(Qtab[0].X, 0);
    mpz_set_ui(Qtab[0].Y, 1);
    mpz_set_ui(Qtab[0].Z, 0);
    Qtab[0].infinity = 1;


    if (!P->infinity) {
        mpz_set(Ptab[1].X, P->X);
        mpz_set(Ptab[1].Y, P->Y);
        mpz_set(Ptab[1].Z, P->Z);
        Ptab[1].infinity = 0;
    } else {Ptab[1].infinity = 1;};

    if (!Q->infinity) {
        mpz_set(Qtab[1].X, Q->X);
        mpz_set(Qtab[1].Y, Q->Y);
        mpz_set(Qtab[1].Z, Q->Z);
        Qtab[1].infinity = 0;
    } else {Qtab[1].infinity = 1;};

    for (unsigned int i = 2; i < M; i++) {
    ec_point_add_proj(&Ptab[i], &Ptab[i - 1], P, E);
    ec_point_add_proj(&Qtab[i], &Qtab[i - 1], Q, E);
}


    // --- Allocation de T[i][j] ---
    ECPointProj **T = malloc(M * sizeof(ECPointProj *));
    for (unsigned int i = 0; i < M; i++) {
        T[i] = malloc(M * sizeof(ECPointProj));
        for (unsigned int j = 0; j < M; j++) {
            ec_point_proj_init(&T[i][j]);
            /*if (Ptab[i].infinity && Qtab[j].infinity) {
                T[i][j].infinity = 1;
                mpz_set_ui(T[i][j].X, 0);
                mpz_set_ui(T[i][j].Y, 1);
                mpz_set_ui(T[i][j].Z, 0);
            } else if (Ptab[i].infinity) {
                mpz_set(T[i][j].X, Qtab[j].X);
                mpz_set(T[i][j].Y, Qtab[j].Y);
                mpz_set(T[i][j].Z, Qtab[j].Z);
                T[i][j].infinity = Qtab[j].infinity;
            } else if (Qtab[j].infinity) {
                mpz_set(T[i][j].X, Ptab[i].X);
                mpz_set(T[i][j].Y, Ptab[i].Y);
                mpz_set(T[i][j].Z, Ptab[i].Z);
                T[i][j].infinity = Ptab[i].infinity;
            } else {*/
                ec_point_add_proj(&T[i][j], &Ptab[i], &Qtab[j], E);
            //}
        }
    }
    
    // Je suis pas sur de l'erreur mais T[0][0] n'a pas la bonne valeur sur mes tests
    //alors que ca devrait valoir 0. Jai essayé de forcé mais ca change pas

    //mpz_set_ui(T[0][0].X, 0);
    //mpz_set_ui(T[0][0].Y, 1);
    //mpz_set_ui(T[0][0].Z, 0);
    //T[0][0].infinity = 1;
    for (unsigned int i = 0; i < M; i++){
            ec_point_proj_clear(&Ptab[i]);
            ec_point_proj_clear(&Qtab[i]);
    }
    free(Ptab);
    free(Qtab);

    return T;
}

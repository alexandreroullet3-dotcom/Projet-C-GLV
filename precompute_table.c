#include "precompute_table.h"
#include "EC_add_proj.h"

/*
 * =========================
 * Pré-calcul de la table pour double multiplication scalaire
 * T[i][j] = i*P + j*Q pour i,j ∈ [0, 2^w - 1]
 * =========================
 */

ECPointProj** precompute_table(const ECPointProj *P, const ECPointProj *Q, 
                            unsigned int w, const ECCurve *E)
{
    unsigned int M = 1U << w; // taille de la fenêtre

    // --- Allocation et init des multiples successifs ---
    ECPointProj *Ptab = malloc(M * sizeof(ECPointProj));
    ECPointProj *Qtab = malloc(M * sizeof(ECPointProj));

    for (unsigned int i = 0; i < M; i++) {
        ec_point_proj_init(&Ptab[i]);
        ec_point_proj_init(&Qtab[i]);
    }

    // Indice 0 = point à l'infini
    Ptab[0].infinity = 1;
    Qtab[0].infinity = 1;

    // Indice 1 = P et Q eux-mêmes si pas infini
    if (!P->infinity) ec_point_proj_copy(&Ptab[1], P); else Ptab[1].infinity = 1;
    if (!Q->infinity) ec_point_proj_copy(&Qtab[1], Q); else Qtab[1].infinity = 1;

    // Calcul des multiples successifs Ptab[i] = i*P, Qtab[j] = j*Q
    for (unsigned int i = 2; i < M; i++) {
        ec_point_add_proj(&Ptab[i], &Ptab[i - 1], P, E); // Ptab[i] = Ptab[i-1] + P
        ec_point_add_proj(&Qtab[i], &Qtab[i - 1], Q, E); // Qtab[i] = Qtab[i-1] + Q
    }

    // --- Allocation et remplissage de la table finale T[i][j] ---
    ECPointProj **T = malloc(M * sizeof(ECPointProj *));
    for (unsigned int i = 0; i < M; i++) {
        T[i] = malloc(M * sizeof(ECPointProj));
        for (unsigned int j = 0; j < M; j++) {
            ec_point_proj_init(&T[i][j]);
            // T[i][j] = Ptab[i] + Qtab[j]
            ec_point_add_proj(&T[i][j], &Ptab[i], &Qtab[j], E);
        }
    }

    // --- Nettoyage des tables temporaires ---
    for (unsigned int i = 0; i < M; i++) {
        ec_point_proj_clear(&Ptab[i]);
        ec_point_proj_clear(&Qtab[i]);
    }
    free(Ptab);
    free(Qtab);

    return T;
}
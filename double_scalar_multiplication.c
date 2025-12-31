#include "double_scalar_multiplication.h"

/*
 * R = kP + lQ
 * Méthode : Shamir avec fenêtre w
 */
void ec_double_scalar_multiplication(
    ECPointProj *R,
    const ECPointProj *P,
    const ECPointProj *Q,
    const mpz_t k,
    const mpz_t l,
    unsigned int w,
    const ECCurve *E)
{
    if (w == 0 || w > 4) return;

    unsigned int M = 1U << w;

    /* =========================
       Pré-calcul des tables
       ========================= */
    
    ECPointProj **T = precompute_table(P, Q, w, E);
    if (!T) return;

    /* =========================
       Initialisation R = O
       ========================= */

    ECPointProj Rtmp;
    ec_point_proj_init(&Rtmp);
    Rtmp.infinity = 1;

    /* =========================
       Taille en bits
       ========================= */

    size_t nb_k = mpz_sizeinbase(k, 2);
    size_t nb_l = mpz_sizeinbase(l, 2);
    size_t n = (nb_k > nb_l) ? nb_k : nb_l;

    size_t d = (n + w - 1) / w;

    /* =========================
       Boucle principale (MSB → LSB)
       ========================= */
    
    for (int i = (int)d - 1; i >= 0; i--) {

    /* R = 2^w * R */
    for (unsigned int j = 0; j < w; j++)
        ec_point_double_proj(&Rtmp, &Rtmp, E);


    unsigned int uk = 0, ul = 0;

    for (unsigned int b = 0; b < w; b++) {
        size_t bit = i * w + b;

        if (bit < nb_k && mpz_tstbit(k, bit))
            uk |= (1U << b);

        if (bit < nb_l && mpz_tstbit(l, bit))
            ul |= (1U << b);
    }

    if (uk || ul)
        ec_point_add_proj(&Rtmp, &Rtmp, &T[uk][ul], E);
}



    /* =========================
       Copie résultat
       ========================= */

    mpz_set(R->X, Rtmp.X);
    mpz_set(R->Y, Rtmp.Y);
    mpz_set(R->Z, Rtmp.Z);
    R->infinity = Rtmp.infinity;


    /* =========================
       Libération mémoire
       ========================= */

    ec_point_proj_clear(&Rtmp);

    for (unsigned int i = 0; i < M; i++) {
        for (unsigned int j = 0; j < M; j++)
            ec_point_proj_clear(&T[i][j]);
        free(T[i]);
    }

    free(T);
}

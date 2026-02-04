#include "double_scalar_multiplication.h"

#include <stdlib.h>

/*
 * =========================
 * Double multiplication scalaire : R = k*P + l*Q
 * Algorithme de Shamir avec fenêtres de taille w
 * =========================
 */
void ec_double_scalar_multiplication(ECPointProj *R,
                                     const ECPointProj *P,
                                     const ECPointProj *Q,
                                     const mpz_t k,
                                     const mpz_t l,
                                     unsigned int w,
                                     const ECCurve *E)
{
    if (w == 0 || w > 4) {
        return; // Failsafe, il ne faut pas w trop grand, sinon le temps de calcul explose
    }

    unsigned int M = 1U << w;

    ECPointProj Pprime, Qprime;
    ec_point_proj_init(&Pprime);
    ec_point_proj_init(&Qprime);
    mpz_t ksgn, lsgn;
    mpz_inits(ksgn, lsgn, NULL);
    
    if (mpz_sgn(k) < 0) {
        ec_point_proj_neg(&Pprime, P);
        mpz_neg(ksgn, k);
    } else {
        ec_point_proj_copy(&Pprime, P);
        mpz_set(ksgn, k);
    }

    if (mpz_sgn(l) < 0) {
        ec_point_proj_neg(&Qprime, Q);
        mpz_neg(lsgn, l);
    } else {
        ec_point_proj_copy(&Qprime, Q);
        mpz_set(lsgn, l);
    }

    /* =========================
       Pré-calcul des tables T[i][j] = i*P + j*Q
       ========================= */
    ECPointProj **T = precompute_table(&Pprime, &Qprime, w, E);
    if (!T) {
        ec_point_proj_clear(&Pprime);
        ec_point_proj_clear(&Qprime);
        mpz_clears(ksgn, lsgn, NULL);
        return;
    }

    /* =========================
       Initialisation Rtmp = O
       ========================= */
    ECPointProj Rtmp;
    ec_point_proj_init(&Rtmp);
    Rtmp.infinity = 1;

    /* =========================
       Taille en bits des scalaires
       ========================= */
    size_t nb_k = mpz_sizeinbase(ksgn, 2);
    size_t nb_l = mpz_sizeinbase(lsgn, 2);
    size_t n = (nb_k > nb_l) ? nb_k : nb_l;
    size_t d = (n + w - 1) / w; // nombre de fenêtres

    /* =========================
       Boucle principale MSB → LSB
       ========================= */
    for (int i = (int)d - 1; i >= 0; i--) {

        // Rtmp = 2^w * Rtmp
        for (unsigned int j = 0; j < w; j++) {
            ec_point_double_proj(&Rtmp, &Rtmp, E);
        }

        // Extraire les w bits correspondants
        unsigned int uk = 0, ul = 0;
        for (unsigned int b = 0; b < w; b++) {
            size_t bit = i * w + b;

            if (bit < nb_k && mpz_tstbit(ksgn, bit)) {
                uk |= (1U << b);
            }

            if (bit < nb_l && mpz_tstbit(lsgn, bit)) {
                ul |= (1U << b);
            }
        }

        ec_point_add_proj(&Rtmp, &Rtmp, &T[uk][ul], E);
    }

    /* =========================
       Copie du résultat final dans R
       ========================= */
    mpz_set(R->X, Rtmp.X);
    mpz_set(R->Y, Rtmp.Y);
    mpz_set(R->Z, Rtmp.Z);
    R->infinity = Rtmp.infinity;

    /* =========================
       Libération mémoire
       ========================= */
    ec_point_proj_clear(&Rtmp);
    ec_point_proj_clear(&Pprime);
    ec_point_proj_clear(&Qprime);
    mpz_clears(ksgn, lsgn, NULL);

    for (unsigned int i = 0; i < M; i++) {
        for (unsigned int j = 0; j < M; j++) {
            ec_point_proj_clear(&T[i][j]);
        }
        free(T[i]);
    }
    free(T);
}

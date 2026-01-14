#include "EC_struct.h"

/*
 * =========================
 * Fonctions pour points affines
 * =========================
 */

/* Initialisation d'un point affine : x et y sont initialisés, point à l'infini par défaut */
void ec_point_affine_init(ECPointAffine *P) {
    mpz_init(P->x);
    mpz_init(P->y);
    P->infinity = 1; /* point à l'infini */
}

/* Libération de la mémoire allouée pour un point affine */
void ec_point_affine_clear(ECPointAffine *P) {
    mpz_clear(P->x);
    mpz_clear(P->y);
}

/* Copie d'un point affine : R <- P */
void ec_point_affine_copy(ECPointAffine *R, const ECPointAffine *P) {
    mpz_set(R->x, P->x);
    mpz_set(R->y, P->y);
    R->infinity = P->infinity;
}

// Compare deux points P et Q.
// Retourne 0 si P == Q (identiques).
// Retourne 1 si P != Q (différents).
int ec_cmp_affine(const ECPointAffine *P, const ECPointAffine *Q) {
    if (P->infinity && Q->infinity) {
        return 0;
    }
    if (P->infinity || Q->infinity) {
        return 1;
    }
    if (mpz_cmp(P->x, Q->x) != 0) {
        return 1;
    }
    if (mpz_cmp(P->y, Q->y) != 0) {
        return 1; 
    }
    return 0;
}

/*
 * =========================
 * Fonctions pour points projectifs
 * =========================
 */

/* Initialisation d'un point projectif : X,Y,Z initialisés, point à l'infini */
void ec_point_proj_init(ECPointProj *P) {
    mpz_init(P->X);
    mpz_init(P->Y);
    mpz_init(P->Z);
    P->infinity = 1; /* point à l'infini */
}

/* Libération de la mémoire allouée pour un point projectif */
void ec_point_proj_clear(ECPointProj *P) {
    mpz_clear(P->X);
    mpz_clear(P->Y);
    mpz_clear(P->Z);
}

/* Copie d'un point projectif : R <- P */
void ec_point_proj_copy(ECPointProj *R, const ECPointProj *P) {
    mpz_set(R->X, P->X);
    mpz_set(R->Y, P->Y);
    mpz_set(R->Z, P->Z);
    R->infinity = P->infinity;
}

/*
 * =========================
 * Fonctions pour la courbe
 * =========================
 */

/* Initialisation des paramètres de la courbe elliptique */
void ec_curve_init(ECCurve *E) {
    mpz_init(E->p);
    mpz_init(E->a);
    mpz_init(E->b);
}

/* Libération de la mémoire associée à la courbe elliptique */
void ec_curve_clear(ECCurve *E) {
    mpz_clear(E->p);
    mpz_clear(E->a);
    mpz_clear(E->b);
}
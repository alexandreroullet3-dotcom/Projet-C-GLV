#include "EC_struct.h"

/*
 * =========================
 * Fonctions pour points affines
 * =========================
 */

/* Initialisation d'un point affine : x et y sont initialisés, point à l'infini par défaut */
void ec_point_affine_init(ECPointAffine *P)
{
    mpz_init(P->x);
    mpz_init(P->y);
    P->infinity = 1; /* point à l'infini */
}

/* Libération de la mémoire allouée pour un point affine */
void ec_point_affine_clear(ECPointAffine *P)
{
    mpz_clear(P->x);
    mpz_clear(P->y);
}

/* Copie d'un point affine : R <- P */
void ec_point_affine_copy(ECPointAffine *R, const ECPointAffine *P)
{
    mpz_set(R->x, P->x);
    mpz_set(R->y, P->y);
    R->infinity = P->infinity;
}


/*
 * =========================
 * Fonctions pour points projectifs
 * =========================
 */

/* Initialisation d'un point projectif : X,Y,Z initialisés, point à l'infini */
void ec_point_proj_init(ECPointProj *P)
{
    mpz_init(P->X);
    mpz_init(P->Y);
    mpz_init(P->Z);
    P->infinity = 1; /* point à l'infini */
}

/* Libération de la mémoire allouée pour un point projectif */
void ec_point_proj_clear(ECPointProj *P)
{
    mpz_clear(P->X);
    mpz_clear(P->Y);
    mpz_clear(P->Z);
}

/* Copie d'un point projectif : R <- P */
void ec_point_proj_copy(ECPointProj *R, const ECPointProj *P)
{
    mpz_set(R->X, P->X);
    mpz_set(R->Y, P->Y);
    mpz_set(R->Z, P->Z);
    R->infinity = P->infinity;
}


/*
 * Conversion affine -> projectif.
 */
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

/*
 * Conversion projectif -> affine.
 */
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


// Compare deux points P et Q.
// Retourne 0 si P == Q (identiques).
// Retourne 1 si P != Q (différents).
int ec_cmp_affine(const ECPointAffine *P, const ECPointAffine *Q)
{
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

int ec_cmp_proj(const ECPointProj *P, const ECPointProj *Q, const ECCurve *E)
{
    if (P->infinity && Q->infinity) {
        return 0;
    }
    if (P->infinity || Q->infinity) {
        return 1;
    }

    mpz_t xp, yp, xq, yq;
    mpz_inits(xp, yp, xq, yq, NULL);
    mpz_mul(xp, P->X, Q->Z);
    mpz_mul(xp, xp, Q->Z);
    mpz_mod(xp, xp, E->p);

    mpz_mul(yp, P->Y, Q->Z);
    mpz_mul(yp, yp, Q->Z);
    mpz_mul(yp, yp, Q->Z);
    mpz_mod(yp, yp, E->p);


    mpz_mul(xq, Q->X, P->Z);
    mpz_mul(xq, xq, P->Z);
    mpz_mod(xq, xq, E->p);
    
    mpz_mul(yq, Q->Y, P->Z);
    mpz_mul(yq, yq, P->Z);
    mpz_mul(yq, yq, P->Z);
    mpz_mod(yq, yq, E->p);
    
    if (mpz_cmp(xp, xq) != 0) {
        mpz_clears(xp, yp, xq, yq, NULL);
        return 1;
    }
    if (mpz_cmp(yp, yq) != 0) {
        mpz_clears(xp, yp, xq, yq, NULL);
        return 1;
    }
    mpz_clears(xp, yp, xq, yq, NULL);
    return 0;
}

/*Tests d'appartenance*/
// Retourne 0 si P n'est pas sur la courbe.
// Retourne 1 si P est sur la courbe.
int is_in_aff(const ECPointAffine *P, const ECCurve *E){
    mpz_t t, s;
    mpz_inits(t, s, NULL);
    mpz_mul(s, P->x, P->x);
    mpz_mul(t, s, P->x);
    mpz_mul(s, s, E->a2);
    mpz_add(t, s, t);
    mpz_mul(s, P->x, E->a);
    mpz_add(t, t, s);
    mpz_add(t, t, E->b);
    mpz_mod(t, t, E->p);
    mpz_mul(s, P->y, P->y);
    mpz_mod(s, s, E->p);
    if (mpz_cmp(s, t)){
        return 1;
    }
    return 0;
}

int is_in_proj(const ECPointProj *P, const ECCurve *E){
    mpz_t t, s;
    mpz_inits(t, s, NULL);

    mpz_mul(s, P->X, P->X);
    mpz_mul(t, s, P->X);

    mpz_mul(s, s, E->a2);
    mpz_mul(s, s, P->Z);
    mpz_mul(s, s, P->Z);
    mpz_add(t, s, t);
    
    mpz_mul(s, P->X, E->a);
    mpz_mul(s, s, P->Z);
    mpz_mul(s, s, P->Z);
    mpz_mul(s, s, P->Z);
    mpz_mul(s, s, P->Z);
    mpz_add(t, t, s);

    mpz_mul(t, P->Z, E->b);
    mpz_mul(s, s, P->Z);
    mpz_mul(s, s, P->Z);
    mpz_mul(s, s, P->Z);
    mpz_mul(s, s, P->Z);
    mpz_mul(s, s, P->Z);

    mpz_mod(t, t, E->p);

    mpz_mul(s, P->Y, P->Y);
    mpz_mul(s, P->Y, P->Y);
    mpz_mod(s, s, E->p);
    if (mpz_cmp(s, t)){
        return 1;
    }
    return 0;
}


/*
 * =========================
 * Fonctions pour la courbe
 * =========================
 */

/* Initialisation des paramètres de la courbe elliptique */
void ec_curve_init(ECCurve *E)
{
    mpz_init(E->p);
    mpz_init(E->a);
    mpz_init(E->b);
    mpz_init_set_ui(E->a2, 0); //Vaut 0 par défaut, non nul seulement dans l'exemple 3
}

/* Libération de la mémoire associée à la courbe elliptique */
void ec_curve_clear(ECCurve *E)
{
    mpz_clear(E->p);
    mpz_clear(E->a);
    mpz_clear(E->b);
    mpz_clear(E->a2);
}
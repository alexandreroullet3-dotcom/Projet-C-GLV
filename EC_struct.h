#ifndef EC_STRUCT_H
#define EC_STRUCT_H

#include <gmp.h>
#include <stdlib.h>


/*
 * =========================
 * Structures de base ECC
 * =========================
 */

/*
 * Point en coordonnées affines :
 *   - (x, y) dans F_p^2
 *   - infinity = 1 représente le point à l'infini O
 */
typedef struct {
    mpz_t x;        /* coordonnée x */
    mpz_t y;        /* coordonnée y */
    int   infinity; /* 1 si point à l'infini, 0 sinon */
} ECPointAffine;

/*
 * Point en coordonnées projectives (Jacobiennes) :
 *   - (X : Y : Z) représente (X/Z^2, Y/Z^3) en affine
 *   - infinity = 1 représente le point à l'infini O
 */
typedef struct {
    mpz_t X;        /* coordonnée X */
    mpz_t Y;        /* coordonnée Y */
    mpz_t Z;        /* coordonnée Z */
    int   infinity; /* 1 si point à l'infini, 0 sinon */
} ECPointProj;

/*
 * Courbe elliptique :
 *   y^2 = x^3 + a x + b (mod p)
 */
typedef struct {
    mpz_t p; /* module premier */
    mpz_t a; /* coefficient a */
    mpz_t b; /* coefficient b */
} ECCurve;

/*
 * =========================
 * Gestion des points affines
 * =========================
 */

/* Initialisation d'un point affine (x,y non définis, point à l'infini) */
void ec_point_affine_init(ECPointAffine *P);

/* Libération de la mémoire associée à un point affine */
void ec_point_affine_clear(ECPointAffine *P);

/* Copie d'un point affine : dst <- src */
void ec_point_affine_copy(ECPointAffine *dst, const ECPointAffine *src);

/*
 * =========================
 * Gestion des points projectifs
 * =========================
 */

/* Initialisation d'un point projectif (X,Y,Z non définis, point à l'infini) */
void ec_point_proj_init(ECPointProj *P);

/* Libération de la mémoire associée à un point projectif */
void ec_point_proj_clear(ECPointProj *P);

/* Copie d'un point projectif : dst <- src */
void ec_point_proj_copy(ECPointProj *dst, const ECPointProj *src);

/*
 * =========================
 * Conversions affine <-> projectif
 * =========================
 */

/* Conversion affine -> projectif (Z = 1 si point non infini) */
void affine_to_proj(ECPointProj *R, const ECPointAffine *P);

/* Conversion projectif -> affine (nécessite une inversion modulo p) */
void proj_to_affine(ECPointAffine *R, const ECPointProj *P, const ECCurve *E);

/*
 * =========================
 * Gestion de la courbe
 * =========================
 */

/* Initialisation des paramètres de la courbe */
void ec_curve_init(ECCurve *E);

/* Libération de la mémoire associée à la courbe */
void ec_curve_clear(ECCurve *E);

#endif /* EC_STRUCT_H */

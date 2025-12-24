#ifndef EC_STRUCT_H
#define EC_STRUCT_H
#include <gmp.h>

// Structure point affine
typedef struct {
    mpz_t x;
    mpz_t y;
    int infinity; //1 si point à l'infini, 0 sinon
} ECPointAffine;

// Structure point projectif
typedef struct {
    mpz_t X;
    mpz_t Y;
    mpz_t Z;
    int infinity; //1 si point à l'infini, 0 sinon
} ECPointProj;


// Structure courbe
typedef struct {
    mpz_t p;
    mpz_t a;
    mpz_t b;
} ECCurve;

// Initialisation / nettoyage
void ec_point_affine_init(ECPointAffine *P);
void ec_point_affine_clear(ECPointAffine *P);
void ec_point_affine_copy(ECPointAffine *dst, const ECPointAffine *src);
void ec_point_proj_init(ECPointProj *P);
void ec_point_proj_clear(ECPointProj *P);
void ec_point_proj_copy(ECPointProj *dst, const ECPointProj *src);
void affine_to_proj(ECPointProj *R, const ECPointAffine *P);
void proj_to_affine(ECPointAffine *R, const ECPointProj *P, const ECCurve *E);
void ec_curve_init(ECCurve *E);
void ec_curve_clear(ECCurve *E);

#endif
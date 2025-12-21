#ifndef EC_H
#define EC_H
#include <gmp.h>

// Structure point affine
typedef struct {
    mpz_t x;
    mpz_t y;
    int infinity; //1 si point à l'infini, 0 sinon
} ECPoint;

// Structure courbe
typedef struct {
    mpz_t p;
    mpz_t a;
    mpz_t b;
} ECCurve;

// Initialisation / nettoyage
void ec_point_init(ECPoint *P);
void ec_point_clear(ECPoint *P);
void ec_curve_init(ECCurve *E);
void ec_curve_clear(ECCurve *E);

#endif
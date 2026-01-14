#ifndef EC_ENDO_PHI_GLV_H
#define EC_ENDO_PHI_GLV_H

#include "EC_square_and_multiply_affine.h"

void ec_endo_phi_point_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t beta);
void ec_endo_phi1_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t beta);
void ec_endo_phi2_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t beta);
void ec_endo_phi3_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t omega);
int mpz_tonnelli_shanks(mpz_t r, const mpz_t n, const mpz_t p) ;
void solve_quadratic_equation(mpz_t r1, mpz_t r2, mpz_t A, mpz_t B, const mpz_t m);
void trouver_constantes_glv(mpz_t beta, mpz_t lambda, const ECCurve *E,
                             const mpz_t n, const ECPointAffine *P, int type);

#endif
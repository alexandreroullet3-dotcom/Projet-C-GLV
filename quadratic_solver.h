#ifndef QUADRATIC_SOLVER_H
#define QUADRATIC_SOLVER_H

#include "EC_struct.h"
#include "EC_square_and_multiply_proj.h"
#include "EC_endo_phi_GLV.h"

// Calcule r tel que r^2 = n mod p
// Retourne 1 si une racine existe, 0 sinon.
int mpz_tonelli_shanks(mpz_t r, const mpz_t n, const mpz_t p);

// Résout l'équation x^2 + A x + B = 0 mod m, retourne deux solutions r1, r2.
void solve_quadratic_equation(mpz_t r1, mpz_t r2, mpz_t A, mpz_t B, const mpz_t m);


// Fonction qui donne beta et lambda (solutions de l'équation caractéristique de phi).
// Prend en entrée la courbe E et le type (1, 2 ou 3).
void trouver_constantes_glv(mpz_t beta, const ECCurve *E, const ECPointProj *P, const mpz_t lambda, int type);


// Fonction utilisée dans le cas 3 pour calculer omega = (1+sqrt(-7))/2.
void find_omega(mpz_t omega, const mpz_t p);

#endif /* QUADRATIC_SOLVER_H */

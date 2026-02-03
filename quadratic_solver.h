#ifndef QUADRA_SOLVER
#define QUADRA_SOLVER
#include "EC_struct.h"
#include "EC_endo_phi_GLV.h"

// Calcule r tel que r^2 = n mod p
// Retourne 1 si une racine existe, 0 sinon.
int mpz_tonelli_shanks(mpz_t r, const mpz_t n, const mpz_t p) ;

// focntion qui Résout une equation x^2 +Ax +B =0 mod m et stock les solution dasn r1 et r2
void solve_quadratic_equation(mpz_t r1, mpz_t r2, mpz_t A, mpz_t B, const mpz_t m);


// Fonction qui nous donne beta et lambda 
// ils sont les solutions de l'equation caracteristique de phi respectrivement modulo p et n
// prend en entrée la courbe E, et le type (1, 2 ou 3)
// cette focntion est a modifier si on rajoute des exemples, car si c'est aucun de c'est trois type,   
// il faut toruver l'équation caracteristique, mais pour cela on doit la connaitre aux préalable
void trouver_constantes_glv(mpz_t beta, const ECCurve *E,
                            int type);


//Fonction utilisée dans le cas 3 pour calculer omega = (1+sqrt(-7))/2
void find_omega(mpz_t omega, const mpz_t p);

#endif
#include "quadratic_solver.h"

/*
 * Tonelli–Shanks
 * Résout r^2 ≡ n (mod p), p premier impair
 * Retourne 1 si solution, 0 sinon
 */
int mpz_tonelli_shanks(mpz_t r, const mpz_t n, const mpz_t p)
{
    /* Cas triviaux */
    if (mpz_cmp_ui(n, 0) == 0) {
        mpz_set_ui(r, 0);
        return 1;
    }

    /* Vérification résidu quadratique */
    if (mpz_legendre(n, p) != 1)
        return 0;
    
    mpz_t tmp;
    mpz_init(tmp);
    /* Cas rapide p ≡ 3 (mod 4) */
    if (mpz_mod_ui(tmp, p, 4) == 3) {
        mpz_t e;
        mpz_init(e);
        mpz_add_ui(e, p, 1);
        mpz_fdiv_q_2exp(e, e, 2);  // (p+1)/4
        mpz_powm(r, n, e, p);
        mpz_clear(e);
        return 1;
    }

    /* Décomposition p-1 = q * 2^s avec q impair */
    mpz_t q;
    mpz_init(q);
    mpz_sub_ui(q, p, 1);

    unsigned long s = 0;
    while (mpz_even_p(q)) {
        mpz_fdiv_q_2exp(q, q, 1);
        s++;
    }

    /* Trouver z non-résidu quadratique */
    mpz_t z;
    mpz_init_set_ui(z, 2);
    while (mpz_legendre(z, p) != -1)
        mpz_add_ui(z, z, 1);

    /* Initialisations */
    mpz_t c, t;
    mpz_inits(c, t, NULL);

    mpz_powm(c, z, q, p);    // c = z^q
    mpz_powm(t, n, q, p);    // t = n^q

    mpz_t e;
    mpz_init(e);
    mpz_add_ui(e, q, 1);
    mpz_fdiv_q_2exp(e, e, 1); // (q+1)/2
    mpz_powm(r, n, e, p);     // r = n^((q+1)/2)

    unsigned long m = s;

    /* Boucle principale */
    while (mpz_cmp_ui(t, 1) != 0) {
        
        mpz_set(tmp, t);

        unsigned long i;
        for (i = 1; i < m; i++) {
            mpz_mul(tmp, tmp, tmp);
            mpz_mod(tmp, tmp, p);
            if (mpz_cmp_ui(tmp, 1) == 0)
                break;
        }

        if (i == m) {
            mpz_clears(q, z, c, t, e, tmp, NULL);
            return 0;
        }

        mpz_t b;
        mpz_init(b);
        mpz_powm_ui(b, c, 1UL << (m - i - 1), p);

        mpz_mul(r, r, b);
        mpz_mod(r, r, p);

        mpz_mul(c, b, b);
        mpz_mod(c, c, p);

        mpz_mul(t, t, c);
        mpz_mod(t, t, p);

        m = i;

        mpz_clear(b);
        mpz_clear(tmp);
    }

    mpz_clears(q, z, c, t, e, NULL);
    return 1;
}

// focntion qui Résout une equation x^2 +Ax +B =0 mod m et stock les solution dasn r1 et r2
void solve_quadratic_equation(mpz_t r1, mpz_t r2, mpz_t A, mpz_t B, const mpz_t m) {
    mpz_t D, sqrt_D, inv_2, t;
    mpz_inits(D, sqrt_D, inv_2, t, NULL);

    // D = A^2 - 4B
    mpz_mul(D, A, A);       
    mpz_set_ui(t, 4);
    mpz_mul(t, t, B);       
    mpz_sub(D, D, t); 
    mpz_mod(D, D, m);

    // Racine carrée de D (utilise ta fonction mpz_sqrtm_prime)
    if (mpz_tonelli_shanks(sqrt_D, D, m) == 0) {
        printf("Erreur : Pas de solution quadratique ! Les entrées étaient\n");
        gmp_printf("A = %Zd, B = %Zd, m = %Zx\n", A, B, m);
        exit(1);
    }

    // Calcul de : (-A ± sqrt) / 2
    mpz_set_ui(inv_2, 2);
    mpz_invert(inv_2, inv_2, m); 
    mpz_neg(t, A);
    
    // r1
    mpz_add(r1, t, sqrt_D);
    mpz_mul(r1, r1, inv_2);
    mpz_mod(r1, r1, m);

    // r2
    mpz_sub(r2, t, sqrt_D);
    mpz_mul(r2, r2, inv_2);
    mpz_mod(r2, r2, m);

    mpz_clears(D, sqrt_D, inv_2, t, NULL);
}

// Fonction qui nous donne beta et lambda 
// ils sont les solutions de l'equation caracteristique de phi respectrivement modulo p et n
// prend en entrée la courbe E, et le type (1, 2 ou 3)
// cette focntion est a modifier si on rajoute des exemples, car si c'est aucun de c'est trois type,   
// il faut toruver l'équation caracteristique, mais pour cela on doit la connaitre aux préalable
void trouver_constantes_glv(mpz_t beta, mpz_t lambda, const ECCurve *E,
                            const mpz_t n, ECPointAffine *P, int type) {
    mpz_t A, B, b1, b2, l1, l2;
    mpz_inits(A, B, b1, b2, l1, l2, NULL);

    // Définissons l'équation , donc les varibale A et B
    if (type == 1) {       // x^2 + x + 1 = 0
        mpz_set_ui(A, 1); 
        mpz_set_ui(B, 1);
    } else if (type == 2) { // x^2 + 1 = 0
        mpz_set_ui(A, 0);
        mpz_set_ui(B, 1);
    } else if (type == 3) { // x^2 - x + 2 = 0
        // Pour type 3, A = -1 mod p
        mpz_sub_ui(A, E->p, 1);
        mpz_set_ui(B, 2);
    }

    // Calculer Beta Modulo p
    solve_quadratic_equation(b1, b2, A, B, E->p);
    
    // On choisit arbitrairement la première solution pour Beta
    mpz_set(beta, b1); 

    if (mpz_cmp_ui(lambda,1)==0){
        printf("passé par la 1\n");
    }
    
    if (P->infinity){
        printf("passé par la 2\n");
    }
    if (mpz_cmp_ui(n, 1) == 0){
        printf("passé par la 3\n");
    }
    
    /*
    // Calculer Lambda (Modulo N) 
    // Ajustement de A pour le modulo N (si type 3, A = n - 1)
    if (type == 3){ 
        mpz_sub_ui(A, n, 1);
    }
    solve_quadratic_equation(l1, l2, A, B, n);
    
    
    // Trouver le bon Lambda (celui lié au beta choisi )
    
    ECPointAffine Q_a, Q_b;
    ec_point_affine_init(&Q_a);
    ec_point_affine_init(&Q_b);

    // On utilise le point P pour tester, car n'importe quelle point de la courbe suffit
    // On calcule Phi(G) avec notre beta choisi
    if (type == 1){
        ec_endo_phi1_affine(&Q_a, P, E, beta);
    }else if (type == 2){
        ec_endo_phi2_affine(&Q_a, P, E, beta);
    }else if (type == 3){ 
        ec_endo_phi3_affine(&Q_a, P, E, beta);
    }
    // Tester avec Lambda 1 : Q = l1 * G
    // ce n'est pas grave de faire d'utiliser square and multiply car ce calcule est fait une seul fois 
    // pour trouver lambda, une fois lambda connu (lambda fait environ la taille de k)
    //on poura utiliser GLV autant de fois que l'on veut sur une même courbe
    ec_scalar_mul_affine(&Q_b, P, l1, E);

    // Comparaison
    if (ec_cmp_affine(&Q_a, &Q_b) == 0) {
        mpz_set(lambda, l1); // C'est l1
    } else {
        mpz_set(lambda, l2); // C'est l2
    }

    // Nettoyage
    mpz_clears(A, B, b1, b2, l1, l2, NULL);
    ec_point_affine_clear(&Q_a);
    ec_point_affine_clear(&Q_b);
    }*/
}

void find_omega(mpz_t omega, const mpz_t p){
    mpz_t a;
    mpz_init(a);
    mpz_set_si(a, -7);
    if (mpz_legendre(a, p) != 1){
        printf("Erreur le nombre premier p ne convient pas, -7 n'est pas résidu quadratique mod p\n");
        return ;
    }
    mpz_tonelli_shanks(a, a, p);
    mpz_add_ui(omega, a, 1);
    mpz_set_ui(a, 2);
    mpz_invert(a, a, p);
    mpz_mul(omega, omega, a);
    mpz_mod(omega, omega, p);
    mpz_clear(a);
}
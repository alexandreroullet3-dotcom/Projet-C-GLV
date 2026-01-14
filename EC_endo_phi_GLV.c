#include "EC_endo_phi_GLV.h"

// fonction qui applique l'endomorphisme phi sur un point P

// fonction qui applique l'endomorphisme phi sur un point P
// il est défini par phi(x,y) = (beta*x mod p, y)  avec beta racine cubiuqe de l'unité
void ec_endo_phi_point_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t beta){
    if(P->infinity){
        Q->infinity =1;
    }else{
    // Q.x =beta *P.x mod p
    mpz_mul(Q->x, P->x, beta);
    mpz_mod(Q->x, Q->x, E->p);

    // Q.y = P.y
    mpz_set(Q->y, P->y);
    Q-> infinity =0;
    }
}

////  out ce qui arrive apres est a tester
// cas 1 :  courbes de la forme y^2 = x^3 + b
// il est défini par phi(x,y) = (beta*x mod p, y)  avec beta racine de l'equation X^2 + X + 1 = 0
void ec_endo_phi1_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t beta){
    if(P->infinity){
        Q->infinity =1;
    }else{
        // Q.x =beta *P.x mod p
        mpz_mul(Q->x, P->x, beta);
        mpz_mod(Q->x, Q->x, E->p);

        // Q.y = P.y
        mpz_set(Q->y, P->y);

        Q-> infinity =0;
    }
}

// cas 2 :  courbes de la forme y^2 = x^3 + ax
// il est défini par phi(x,y) = ( -x , beta*y)  avec beta racine de l'equation X^2 + 1 = 0
void ec_endo_phi2_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t beta){
    if(P->infinity){
        Q->infinity =1;
    }else{
        // Q.x = -P.x mod p
        mpz_sub(Q->x, E->p, P->x);

        // Q.y = beta * P.y mod p
        mpz_mul(Q->y, P->y, beta);
        mpz_mod(Q->y, Q->y, E->p);

        Q-> infinity =0;
    }
}

// cas 3 :  courbes de la forme y^2 = x^3 + -3/4 x^2 - 2x - 1
// omega est defini par (1 + sqrt(-7))/2 et a par (omega - 3)/4 (elle seront calculer dasn Fp au préalable)
// il est défini par 
//phi(x,y) = (omega^(-2) * (x^2 - omega))/(x - a) , omega^(-3) * y * (x^2 - 2ax + omega))/(x - a)^2 ) 
void ec_endo_phi3_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t omega){
    if(P->infinity){
        Q->infinity =1;
    }else{
        mpz_t a, inv, inv_w2, inv_w3, num, inv_den, x2, inv_4;
        mpz_inits(a, inv, inv_w2, inv_w3, num, inv_den, x2, inv_4, NULL);

        mpz_set_ui(inv_4, 4);
        mpz_invert(inv_4, inv_4, E->p);
        mpz_sub_ui(a, omega, 3); // a = omega - 3
        mpz_mul(a, a, inv_4); // a = (omega-3) * 4^-1
        mpz_mod(a, a, E->p); // a mod p

        // Calcul de inv= 1 / omega
        mpz_invert(inv, omega, E->p);

        // Calcul de inv_w2 = 1 / (omega^2)
        mpz_mul(inv_w2, inv, inv);  
        mpz_mod(inv_w2, inv_w2, E->p);

        // Calcul de inv_w3 = 1 / (omega^3)
        mpz_mul(inv_w3, inv_w2, inv); 
        mpz_mod(inv_w3, inv_w3, E->p); 

        // Calcul de 1 / (x - a)
        mpz_sub(inv_den, P->x, a);
        mpz_mod(inv_den, inv_den, E->p);

        if (mpz_invert(inv_den, inv_den, E->p) == 0) {
        Q->infinity = 1; // C'est ici qu'on capture l'erreur mathématique
        } 
        else {
            // Calcul x^2
            mpz_mul(x2, P->x, P->x);
            mpz_mod(x2, x2, E->p);

            // Q.x = inv_w2 * (P.x^2 - omega) * inv_den
            mpz_sub(num, x2, omega);
            mpz_mul(Q->x, num, inv_den);
            mpz_mul(Q->x, Q->x, inv_w2);
            mpz_mod(Q->x, Q->x, E->p);

            // Q.y = inv_w3 * P.y * (Px^2 - 2aP.x + omega) * inv_den^2
            mpz_mul(Q->y, a, P->x);
            mpz_mul_ui(Q->y, Q->y, 2);
            mpz_sub(num, x2, Q->y);
            mpz_add(num, num, omega);
            mpz_mul(inv_den, inv_den, inv_den); // (1/(x-a))^2
            mpz_mod(inv_den, inv_den, E->p);
            mpz_mul(Q->y, num, inv_den);      
            mpz_mul(Q->y, Q->y, P->y);  
            mpz_mul(Q->y, Q->y, inv_w3);    
            mpz_mod(Q->y, Q->y, E->p);
            Q-> infinity =0;   
        }
        mpz_clears(a, inv, inv_w2, inv_w3, num, inv_den, x2, inv_4, NULL); 
    }
}

// Calcule r tel que r^2 = n mod p
// Retourne 1 si une racine existe, 0 sinon.
int mpz_tonnelli_shanks(mpz_t r, const mpz_t n, const mpz_t p) {
    // Cas triviaux
    if (mpz_cmp_ui(n, 0) == 0) {
        mpz_set_ui(r, 0);
        return 1;
    }
    
    // Vérification Résidu Quadratique
    if (mpz_legendre(n, p) != 1) {
        return 0; 
    }

    // Cas simple : p = 3 mod 4
    if (mpz_tstbit(p, 0) == 1 && mpz_tstbit(p, 1) == 1){ 
        mpz_t exp;
        mpz_init(exp);
        mpz_add_ui(exp, p, 1);
        mpz_divexact_ui(exp, exp, 4); 
        mpz_powm(r, n, exp, p); 
        mpz_clear(exp);
        return 1;
    }

    // Algorithme de Tonelli-Shanks pour p = 1 mod 4
    mpz_t s, q, z, c, t, m, b, tp, tp2, i, temp_t;
    mpz_inits(s, q, z, c, t, m, b, tp, tp2, i, temp_t, NULL);

    // Écrire p - 1 = q * 2^s
    mpz_sub_ui(q, p, 1);
    mpz_set_ui(s, 0);
    while (mpz_even_p(q)) {
        mpz_divexact_ui(q, q, 2);
        mpz_add_ui(s, s, 1);
    }

    // Trouver z un non-résidu de maniere aléatoire
    mpz_set_ui(z, 2);
    while (mpz_legendre(z, p) != -1) {
        mpz_add_ui(z, z, 1);
    }

    // Initialisation
    mpz_powm(c, z, q, p);
    
    mpz_add_ui(r, q, 1);
    mpz_divexact_ui(r, r, 2); 
    mpz_powm(r, n, r, p); // r initial

    mpz_powm(t, n, q, p); // t initial
    mpz_set(m, s); // m initial

    // Boucle principale
    while (mpz_cmp_ui(t, 1) != 0) {
        // Trouver i tel que t^(2^i) = 1
        mpz_set_ui(i, 1);
        
        int found = 0;
        while (mpz_cmp(i, m) < 0) {
            unsigned long i_ul = mpz_get_ui(i);
            
            // Calcul de t^(2^i) avec mpz_powm
            mpz_set_ui(tp, 1); 
            mpz_mul_2exp(tp, tp, i_ul); // tp = 2^i
            mpz_powm(temp_t, t, tp, p); // temp_t = t^tp mod p
            
            if (mpz_cmp_ui(temp_t, 1) == 0) {
                found = 1;
                break;
            }
            mpz_add_ui(i, i, 1);
        }

        if (!found) {
            mpz_clears(s, q, z, c, t, m, b, tp, tp2, i, temp_t, NULL);
            return 0;
        }

        // b = c^(2^(m-i-1))
        mpz_sub(tp, m, i);
        mpz_sub_ui(tp, tp, 1); // tp = m - i - 1
        
        unsigned long exp_b = mpz_get_ui(tp);
        
        // Calcul de l'exposant 2^exp_b
        mpz_set_ui(tp2, 1);
        mpz_mul_2exp(tp2, tp2, exp_b); // tp2 = 2^(m-i-1)
        
        // Calcul final de b
        mpz_powm(b, c, tp2, p); // b = c^tp2 mod p

        // Mise à jour r, c, t, m
        mpz_mul(r, r, b);
        mpz_mod(r, r, p); 

        mpz_mul(c, b, b);
        mpz_mod(c, c, p); 

        mpz_mul(t, t, c);
        mpz_mod(t, t, p); 

        mpz_set(m, i); 
    }

    mpz_clears(s, q, z, c, t, m, b, tp, tp2, i, temp_t, NULL);
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
    if (mpz_tonnelli_shanks(sqrt_D, D, m) == 0) {
        printf("Erreur : Pas de solution quadratique !\n");
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
                            const mpz_t n, const ECPointAffine *P, int type) {
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
}
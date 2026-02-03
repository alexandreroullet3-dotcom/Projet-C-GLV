#include "EC_endo_phi_GLV.h"

// fonction qui applique l'endomorphisme phi sur un point P


// cas 1 :  courbes de la forme y^2 = x^3 + b
// l'endomorphisme est défini par phi(x,y) = (beta*x mod p, y)  avec beta racine de l'equation X^2 + X + 1 = 0
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
// l'endomorphisme est défini par phi(x,y) = ( -x , beta*y)  avec beta racine de l'equation X^2 + 1 = 0
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
// beta est defini par (1 + sqrt(-7))/2 et a par (beta - 3)/4 (elle seront calculer dasn Fp au préalable)
// l'endomorphisme est défini par 
//phi(x,y) = (beta^(-2) * (x^2 - beta))/(x - a) , beta^(-3) * y * (x^2 - 2ax + beta))/(x - a)^2 ) 
void ec_endo_phi3_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t beta){
    if(P->infinity){
        Q->infinity =1;
    }else{
        mpz_t a, inv, inv_b2, inv_b3, num, inv_den, x2, inv_4;
        mpz_inits(a, inv, inv_b2, inv_b3, num, inv_den, x2, inv_4, NULL);

        mpz_set_ui(inv_4, 4);
        mpz_invert(inv_4, inv_4, E->p);
        mpz_sub_ui(a, beta, 3); // a = beta - 3
        mpz_mul(a, a, inv_4); // a = (beta-3) * 4^-1
        mpz_mod(a, a, E->p); // a mod p

        // Calcul de inv= 1 / beta
        mpz_invert(inv, beta, E->p);

        // Calcul de inv_b2 = 1 / (beta^2)
        mpz_mul(inv_b2, inv, inv);  
        mpz_mod(inv_b2, inv_b2, E->p);

        // Calcul de inv_b3 = 1 / (beta^3)
        mpz_mul(inv_b3, inv_b2, inv); 
        mpz_mod(inv_b3, inv_b3, E->p); 

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

            // Q.x = inv_b2 * (P.x^2 - beta) * inv_den
            mpz_sub(num, x2, beta);
            mpz_mul(Q->x, num, inv_den);
            mpz_mul(Q->x, Q->x, inv_b2);
            mpz_mod(Q->x, Q->x, E->p);

            // Q.y = inv_b3 * P.y * (Px^2 - 2aP.x + beta) * inv_den^2
            mpz_mul(Q->y, a, P->x);
            mpz_mul_ui(Q->y, Q->y, 2);
            mpz_sub(num, x2, Q->y);
            mpz_add(num, num, beta);
            mpz_mul(inv_den, inv_den, inv_den); // (1/(x-a))^2
            mpz_mod(inv_den, inv_den, E->p);
            mpz_mul(Q->y, num, inv_den);      
            mpz_mul(Q->y, Q->y, P->y);  
            mpz_mul(Q->y, Q->y, inv_b3);    
            mpz_mod(Q->y, Q->y, E->p);
            Q-> infinity =0;   
        }
        mpz_clears(a, inv, inv_b2, inv_b3, num, inv_den, x2, inv_4, NULL); 
    }
}
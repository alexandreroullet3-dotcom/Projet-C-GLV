#include "EC_endo_phi_GLV.h"


// Division avec arrondi à l'entier le plus proche nécessaire dans le calcule de 1 et b2
// result = round(numeratuer  / denominateur)
// je met static devant car elle ne sera jamais réutiliser en dehors d'ici
static void mpz_div_arr_proch(mpz_t result, const mpz_t num, const mpz_t den){
    mpz_t r, q, abs_den, seuil;
    mpz_inits(r, q, abs_den, seuil, NULL);

    // Div euclidienne num = q * den + r
    mpz_fdiv_qr(q, r, num, den);

    // seuil = |den| / 2
    mpz_abs(abs_den, den);
    mpz_div_ui(seuil, abs_den, 2); 
    
    // on prend la valeur absolue du reste
    mpz_abs(r, r); 

    // Si le reste est supérieur à la moitié du diviseur, on doit modidfier le quotient
    if (mpz_cmp(r, seuil) > 0) {
        // Si num et den sont de même signe, le quotient est positif donc on ajoute 1
        if ((mpz_sgn(num) > 0) == (mpz_sgn(den) > 0)) {
            mpz_add_ui(q, q, 1);
        // Sinon, le quotient est négatif on soustrait 1    
        } else {
            mpz_sub_ui(q, q, 1);
        }
    }
    mpz_set(result, q);
    mpz_clears(r, q, abs_den, seuil, NULL);
}


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


// decompose k en k1 et k2 tel que k = (k1 + k2 * lambda) mod n avec lambda valeur propre de l'endomophisme phi
// x1, x2,y1,y2 les vecteur de base du réseau, 
// calculer à l'aide de l'algorithme d'euclide apliquer sur n et lambda
// v1=(x1,y1) et v2=(x2,y2) et v=(k,0)
void ec_glv_decompose(mpz_t k1, mpz_t k2, const mpz_t k, const mpz_t x1, const mpz_t y1, 
                    const mpz_t x2, const mpz_t y2){
    mpz_t det, b1, b2, num, t;
    mpz_inits(det, b1, b2, num, t, NULL);

    // résoudre le système linéaire k= b1*x1 + b2*x2 et 0= b1*y1 + b2*y2
    // On veut approcher le vecteur (k, 0) par b1*v1 + b2*v2
    // det = x1*y2 - x2*y1
    mpz_mul(det, x1, y2);
    mpz_mul(t, x2, y1);
    mpz_sub(det, det, t);

    // on utlise la methode de Cramer pour résoudre le système
    // b1 = round( (k * y2) / det )
    mpz_mul(num, k, y2);
    mpz_div_arr_proch(b1, num, det);

    // b2 = round( (-k * y1) / det )
    mpz_mul(num, k, y1);
    mpz_neg(num, num);
    mpz_div_arr_proch(b2, num, det);

    // k1 = k - b1*x1 - b2*x2
    mpz_mul(t, b1, x1);
    mpz_sub(k1, k, t);  
    mpz_mul(t, b2, x2);
    mpz_sub(k1, k1, t); 

    // k2 = -b1*y1 - b2*y2
    mpz_mul(k2, b1, y1);
    mpz_neg(k2, k2);       
    mpz_mul(t, b2, y2);
    mpz_sub(k2, k2, t); 

    // Nettoyage
    mpz_clears(det, b1, b2, num, t, NULL);
}
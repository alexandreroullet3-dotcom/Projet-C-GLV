#include "glv_decompose.h"

void mpz_round_div(mpz_t result, const mpz_t num, const mpz_t den)
{
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

void glv_nearest_vector(Z2 *v,
                        const mpz_t k,
                        const Z2 *v1,
                        const Z2 *v2, const mpz_t n)
{   mpz_t x1, x2, y1, y2, k1, k2;
    mpz_inits(x1, x2, y1, y2, k1, k2, NULL);
    mpz_set(x1, v1->x);
    mpz_set(y1, v1->y);
    mpz_set(x2, v2->x);
    mpz_set(y2, v2->y);
    mpz_t det, b1, b2, num, t;
    mpz_inits(det, b1, b2, num, t, NULL);

    // résoudre le système linéaire k = b1*x1 + b2*x2 et 0 = b1*y1 + b2*y2
    // On veut approcher le vecteur (k, 0) par b1*v1 + b2*v2
    // det = x1*y2 - x2*y1
    mpz_mul(det, x1, y2);
    mpz_mul(t, x2, y1);
    mpz_sub(det, det, t);

    // on utlise la methode de Cramer pour résoudre le système
    // b1 = round( (k * y2) / det )
    mpz_mul(num, k, y2);
    mpz_round_div(b1, num, det);

    // b2 = round( (-k * y1) / det )
    mpz_mul(num, k, y1);
    mpz_neg(num, num);
    mpz_round_div(b2, num, det);

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
    

    mpz_mod(k1, k1, n);
    mpz_mod(k2, k2, n);
    mpz_set(v->x,k1);
    mpz_set(v->y,k2);
    
    mpz_t lambda;
    mpz_init_set_str(lambda, "81e16b4d3131f1322cf0ab2ba439286b3962df578bf08b5dc04a37e1c8e31333", 16);

    mpz_t tmp;
    mpz_init(tmp);

    mpz_mul(tmp, v->y, lambda);
    mpz_add(tmp, tmp, v->x);
    mpz_mod(tmp, tmp, n);

    if (mpz_cmp(tmp, k) != 0)
        //{printf("DECOMP FAUSSE\n");}

    // Nettoyage
    mpz_clears(det, b1, b2, k1, k2, x1, x2, y1, y2, num, t, lambda, tmp, NULL);
    
}
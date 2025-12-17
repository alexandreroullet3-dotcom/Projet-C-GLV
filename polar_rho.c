#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

void c_alea(mpz_t c, mpz_t n, gmp_randstate_t etat) {
    mpz_urandomm(c, etat, n);
}

void function_f(mpz_t out, mpz_t in, mpz_t c, mpz_t n){
    mpz_mul(out, in, in);
    mpz_add(out, out, c);
    mpz_mod(out, out, n);
}

void pollar_rho(mpz_t res, mpz_t n){
    mpz_t c, x, y, d, diff;
    gmp_randstate_t etat;
    mpz_inits(c, x, y, d, diff, NULL);
    gmp_randinit_default(etat);
    gmp_randseed_ui(etat, time(NULL));
    while(1){
        mpz_set_ui(x, 2);
        mpz_set_ui(y, 2);
        mpz_set_ui(d, 1);
        c_alea(c, n, etat);
        while(mpz_cmp_ui(d, 1) == 0){
            function_f(x, x, c, n);

            function_f(y, y, c, n);
            function_f(y, y, c, n);

            mpz_sub(diff, x, y);
            mpz_abs(diff, diff);
            mpz_gcd(d, diff, n);
        }

        if(mpz_cmp(d, n) != 0){
            mpz_set(res, d);
            break;
        }
    }
    mpz_clears(c, x, y, d, diff, NULL);
    gmp_randclear(etat);
}

int main(){
    mpz_t n, res;
    mpz_inits(n, res, NULL);
    mpz_set_str(n, "2177241218019392284455749961185783753335013327591", 10); 
    pollar_rho(res, n);
    gmp_printf("Un facteur trouve : %Zd\n", res);
    mpz_clears(n, res, NULL);
    return 0;
}
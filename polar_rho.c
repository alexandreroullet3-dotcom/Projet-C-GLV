#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

void c_alea(mpz_t c, mpz_t n, gmp_randstate_t etat) {
    mpz_urandomm(c, etat, n);
}

void function_f(mpz_t y,mpz_t x,mpz_t c,mpz_t n){
    mpz_pow_ui(y,x,2);
    mpz_add (y,y,c);
    mpz_mod(y, x, n)
}

void polar_rho(mpz_t res,mpz_t n){
    mpz_t c,x,y,d,abs;
    gmp_randstate_t etat;
    mpz_inits(c,x,y,d,abs);
    while(1){
        mpz_set_ui(x, 2);
        mpz_set_ui(y, 2);
        mpz_set_ui(d, 1);
        gmp_randinit_default(etat);
        gmp_randseed_ui(etat, time(NULL));
        alea(c,n,etat);
        while(mpz_cmp_ui(d,1)==0){
            function_f(x,x,c,n);
            function_f(y,x,c,n);
            gmp_sub(abs,y,x);
            mpz_gcd(d,abs);
        }
        if(mpz_cmp(d,n)!=0){
            mpz_set(res,d);
            break;
        }
    }
    mpz_clears(c,n,x,y,d,abs,NULL);
    gmp_randclear(state);
}

int main(){
    mpz_t n, res;
    mpz_inits(n, f, NULL);
    mpz_set_str(n, "199", 10);
    pollard_rho(res, n);
    gmp_printf("Un facteur trouve : %Zd\n", f);
    mpz_clears(n, f, NULL);
    return 1;
}
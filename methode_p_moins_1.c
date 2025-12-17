#include <stdio.h>



int sq_mult(int x, unsigned int d){
    int temp=1; int puiss=x;
    while (d!=0){
        if (d&1){
            temp = temp*puiss;
        }
        puiss = puiss*puiss;
        d=d>>1;
    } 
    return temp;
}

int algorithme_d_euclide(int a, int b){
    int r;
    r = a%b;
    if (r==0){
        return b;
    };
    return algorithme_d_euclide(b, r);
}

int factoriser_n (int n){
    int a = 2;
    for (int i = 2; i < n; i++){
        int r = algorithme_d_euclide(a-1, n);
        if ((r != 1) && (r != n)){
            return r;
        }
        a = sq_mult(a, i);
    }
    return 1;
}

int main(){
    int n = 299;
    printf("un facteur non trivial est %d\n", factoriser_n(n));
    return 0;
}
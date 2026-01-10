#include <stdlib.h>
#include <gmp.h>
#include "EC_struct.h"
#include "EC_add_proj.h"
#include "EC_endo_phi_GLV.h"
#include "precompute_table.h"
#include "EC_conversions.h"

// Multiplication scalaire GLV  R = k * P

void ec_scal_mul_glv(ECPointProj *R, const ECPointProj *P, const mpz_t k, const ECCurve * E, const mpz_t x1, 
                    const mpz_t y1, const mpz_t x2, const mpz_t y2, const mpz_t beta){
    // on commence par decomposer k
    mpz_t k1,k2;
    mpz_inits(k1,k2,NULL);
    ec_glv_decompose(k1,k2,k,x1,y1,x2,y2);

    // on va creer les tables de précalcule
    ECPointProj P2,Q2;
    ec_point_proj_init(&P2);
    ec_point_proj_init(&Q2);

    ec_point_proj_copy(&P2, P);// on copie juste P dans P2
    //Gestion du signe , si K1<0 alors on prend -P pour traviller (donc P2 = -P)
    if(mpz_sgn(k1) < 0){ // (X,-Y,Z)
        mpz_abs(k1,k1);
        mpz_sub(P2.Y,E->p,P2.Y);
    }

    // Calculons Q = phi(P) 
    //pour cela on passe en affine pour simplifier les calcules et appliquer seulement x*beta
    ECPointAffine P_aff, Q_aff;
    ec_point_affine_init(&P_aff);
    ec_point_affine_init(&Q_aff);
    proj_to_affine(&P_aff, P,E);
    ec_endo_phi_point_affine(&Q_aff,&P_aff,E,beta);

    //on regarde le signe de k2, si il est négatif on prend -Q
    if(mpz_sgn(k2) < 0){
        mpz_abs(k2,k2);
        mpz_sub(Q_aff.y,E->p,Q_aff.y);
    }

    // on repase en projectif pour simplifier les calcule et eviter les division
    affine_to_proj(&Q2,&Q_aff);

    // passons auprécalcul de la table
    int m=2;
    ECPointProj **T =precompute_table( &P2, &Q2, m, E);

    // on regarde combien de bit il faut pour ecrire k1 et k2 en binaire et on prend le plus grand des deux
    size_t nb_bit_k1 = mpz_sizeinbase(k1,2);
    size_t nb_bit_k2 = mpz_sizeinbase(k2,2);
    int max_bit =(nb_bit_k1 > nb_bit_k2) ? nb_bit_k1 : nb_bit_k2;

    // on initialise R à l'infini
    R->infinity = 1;

    // lisse les tailles 
    if(max_bit % m != 0) {
        max_bit += ( m - (max_bit % m));
    }

    for(int i = max_bit - 1; i >= 0; i -= m){
        //calcul R = 2^m * R
        for(int j = 0; j < m; j++){
            ec_point_double_proj(R,R,E);
        }

        int u = 0;
        int v = 0;

        for(int h = 0; h < m; h++){
            int index = i - h;
            if ( index >= 0){
                // regarde le bit numero index de k1 et renvoie 1 si il vaut 1 t 0 sinon
                if( mpz_tstbit(k1, index)) u |= (1<< (m - 1 - h));
                // pareil pour k2
                if( mpz_tstbit(k2, index)) v |= (1<< (m - 1 - h));
            }
        }

        // on ajoute les point précalculé T[u][v]
        // si u=0 et v=0, on ajoute l'infini, donc ca ne fait rien
        ec_point_add_proj(R, R, &T[u][v], E);
    }    

    // on nettoie toute la memoire
    int M = 1<< m;
    for( int l =0 ; l<M; l++){
        for( int f =0 ; f<M; f++){
            ec_point_proj_clear(&T[l][f]);
        }
        free(T[l]);
    }
    free(T);
    mpz_clears(k1, k2, NULL );
    ec_point_proj_clear(&P2);
    ec_point_proj_clear(&Q2);
    ec_point_affine_clear(&P_aff);
    ec_point_affine_clear(&Q_aff);
}

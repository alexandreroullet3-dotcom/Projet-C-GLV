#ifndef EC_ENDO_PHI_GLV_H
#define EC_ENDO_PHI_GLV_H

#include "EC_struct.h"
/*
 * Endomorphismes GLV pour différents types de courbes.
 */
void ec_endo_phi1_affine(ECPointAffine *Q,
                         const ECPointAffine *P,
                         const ECCurve *E,
                         const mpz_t beta);
void ec_endo_phi2_affine(ECPointAffine *Q,
                         const ECPointAffine *P,
                         const ECCurve *E,
                         const mpz_t beta);
void ec_endo_phi3_affine(ECPointAffine *Q,
                         const ECPointAffine *P,
                         const ECCurve *E,
                         const mpz_t omega);

#endif
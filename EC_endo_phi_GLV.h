#ifndef EC_ENDO_PHI_GLV_H
#define EC_ENDO_PHI_GLV_H

#include "EC_struct.h"

void ec_endo_phi_point_affine(ECPointAffine *Q, const ECPointAffine *P, const ECCurve *E, const mpz_t beta);

#endif
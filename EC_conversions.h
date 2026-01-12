#include "EC_struct.h"

void affine_to_proj(ECPointProj *R, const ECPointAffine *P);

void proj_to_affine(ECPointAffine *R,
                    const ECPointProj *P,
                    const ECCurve *E);

/*
 * =========================
 * Comparaisons affine et projectif
 * =========================
 */
int ec_cmp_aff(ECPointAffine *P, ECPointAffine *Q);

int ec_cmp_proj(ECPointProj *P, ECPointProj *Q, ECCurve *E);
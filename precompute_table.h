#ifndef PRECOMPUTE_TABLE_H
#define PRECOMPUTE_TABLE_H

#include "EC_struct.h"

/*
 * =========================
 * Pré-calcul pour double multiplication scalaire
 * =========================
 *
 * T[i][j] = i*P + j*Q pour i,j ∈ [0, 2^w - 1]
 * w = taille de la fenêtre (1..4)
 *
 * Retourne un tableau 2D alloué dynamiquement :
 *   - Chaque ECPointProj doit être libéré avec ec_point_proj_clear avant free
 */
ECPointProj **precompute_table(const ECPointProj *P, const ECPointProj *Q,
                               unsigned int w, const ECCurve *E);

#endif /* PRECOMPUTE_TABLE_H */

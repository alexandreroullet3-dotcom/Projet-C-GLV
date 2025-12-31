#ifndef PRECOMPUTE_TABLE_H
#define PRECOMPUTE_TABLE_H

#include "EC_struct.h"

/**
 * @brief Pré-calculer la table T[i][j] = i*P + j*Q pour double scalar multiplication
 *
 * @param P Le premier point
 * @param Q Le second point
 * @param w La taille de la fenêtre (1..4)
 * @param E La courbe elliptique
 * @return ECPointProj** Tableau 2D alloué dynamiquement de taille 2^w x 2^w
 *         Chaque point doit être libéré avec ec_point_proj_clear avant free.
 */
ECPointProj **precompute_table(const ECPointProj *P,
                                      const ECPointProj *Q,
                                      unsigned int w,
                                      const ECCurve *E);

#endif

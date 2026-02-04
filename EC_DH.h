#ifndef EC_DH_H
#define EC_DH_H

#include "EC_GLV.h"
#include "glv_curves.h"

/*
 * Génération de clés Diffie-Hellman avec GLV.
 */
void ec_dh(GLVCurve *curve);

#endif

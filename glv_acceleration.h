#ifndef GLV_ACCELERATION_H
#define GLV_ACCELERATION_H

#include "EC_GLV.h"
#include "glv_curves.h"

/*
 * Mesure l'accélération GLV sur une courbe donnée.
 */
void glv_acceleration(const GLVCurve *curve, int N_tests);

#endif
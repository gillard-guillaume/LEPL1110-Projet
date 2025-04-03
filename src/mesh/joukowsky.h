#ifndef JOUKOWSKY_H
#define JOUKOWSKY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
# include "../fem/fem.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double complex joukowskyTransform(double complex zeta);

double complex inverseJoukowskyTransform(double complex z);

int joukowsky(femGeo *theGeometry);

#endif
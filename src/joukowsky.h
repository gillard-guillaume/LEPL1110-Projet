#ifndef JOUKOWSKY_H
#define JOUKOWSKY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double complex joukowskyTransform(double complex zeta);

int joukowsky(double R, double mu_x, double mu_y, int N);

#endif

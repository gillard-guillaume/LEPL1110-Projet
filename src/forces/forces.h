#ifndef FORCES_H
#define FORCES_H

#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double f(double x, double L, double b);
double g(double x, double L, double b);
double foilProfile(double x, double y, double lenght, double lift);

#endif // FORCES_H
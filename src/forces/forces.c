#include "forces.h"
#include <math.h> 

double f(double x, double L, double b) {
    double xi = 2.0 * x / b;
    double arg = 1.0 - xi * xi;
    if (arg < 0.0) arg = 0.0; // Security
    return (L / b) * 4.0 * M_PI * sqrt(arg);
}

double g(double x, double L, double b) {
    double xi = 2.0 * x / b;
    if (xi > 1.0) xi = 1; // Security
    return f(0.0, L, b) * (1.0 - xi);
}

double foilProfile(double x, double y, double length, double lift) {
    double L1 = lift;   // Amplitude haut
    double L2 = lift;  // Amplitude bas
    double b = length / 2.0;

    if (y > 0)
        return (x <= 0) ? -f(x, L1, b) : -g(x, L1, b);  // Extrados
    else
        return (x <= 0) ? f(x, L2, b) : g(x, L2, b);  // Intrados
}
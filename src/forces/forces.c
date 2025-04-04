#include "forces.h"
#include <math.h>  // Pour sqrt

double f(double x, double L, double b) {
    // fct "portance racine" sur l'extrados/intrados avant
    double xi = 2.0 * x / b;
    double arg = 1.0 - xi * xi;
    if (arg < 0.0) arg = 0.0; // Sécurité pour éviter sqrt négatif
    return (L / b) * 4.0 * M_PI * sqrt(arg);
}

double g(double x, double L, double b) {
    // fct "portance linéaire" sur extrados/intrados arrière
    // Correction ici pour ne PAS diviser par zéro
    if (fabs(x) < 1e-8) x = (x >= 0) ? 1e-8 : -1e-8;  // éviter 1/0
    return f(0.0, L, b) * (1.0 - 2.0 * x / b);
}

double foilProfile(double x, double y, double length, double lift) {
    double L1 = lift;   // Amplitude haut
    double L2 = lift;  // Amplitude bas
    double b = length / 2.0;

    if (y > 0)
        return (x <= 0) ? f(x, L1, b) : g(x, L1, b);  // Extrados
    else
        return (x <= 0) ? f(x, L2, b) : g(x, L2, b);  // Intrados
}
#ifndef JOUKOWSKY_H
#define JOUKOWSKY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// Function to apply the Joukowsky transformation
double complex joukowsky(double complex zeta);

// Function to generate airfoil points and save to a file
void generate_airfoil(double R, double mu_x, double mu_y, int N, const char *filename);

// Function to generate circles at three locations on the airfoil
void generate_circles(int N, const char *airfoil_file, const char *circle1_file, 
                      const char *circle2_file, const char *circle3_file);

#endif // JOUKOWSKY_H

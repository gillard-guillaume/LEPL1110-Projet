# ifndef MESH_H
# define MESH_H

# include <stdio.h>
# include <stdlib.h>
#include "../gmsh/gmsh-4.13.1-Linux64-sdk/include/gmshc.h"


double geoSize(double x, double y);

int generateSurface(double *x, double *y, int N);

int getPoints (const char *filename, double **x, double **y, int *N);

int wing();


#endif
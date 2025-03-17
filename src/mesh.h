# ifndef MESH_H
# define MESH_H

# include <stdio.h>
# include <stdlib.h>
# include <complex.h>
# include "fem.h"
#include "joukowsky.h"
#include "../gmsh/gmsh-4.13.1-Linux64-sdk/include/gmshc.h"


double geoSize(double x, double y);

int generateSurface(double *x, double *y, int N);

int getPoints(double *x, double *y, int n, char *filename);

int wing(femGeo *theGeometry);


#endif
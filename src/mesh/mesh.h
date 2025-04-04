# ifndef MESH_H
# define MESH_H

# include <stdio.h>
# include <stdlib.h>
# include <complex.h>
# include "../fem/fem.h"
#include "joukowsky.h"
#include "gmshc.h"

femGeo *geoGetGeometry();

double geoSizeDefault(double x, double y);

double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data);

double interpolationH(double x, double x0, double x1, double h0, double h);

double geoSize(double x, double y);

void geoMeshGenerate();

void geoInitialize();

void geoFinalize();

void geoSetSizeCallback(double (*geoSize)(double x, double y));

void geoSetDomainName(int iDomain, char *name);

int geoGetDomain(char *name);

void geoMeshImport();

void geoMeshPrint();

void geoMeshWrite(const char *filename);

void geoMeshRead(const char *filename);

void geoMeshUnfuck();


#endif
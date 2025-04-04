#include "generateMesh.h"

void generateMesh(double R, double mu_x, double mu_y){
    clock_t start = clock();
    int N = 100;
    double h = 0.1;

    femGeo* theGeometry = geoGetGeometry();    
    theGeometry->h            =  0.1;  
    theGeometry->R            =  R;
    theGeometry->muX         =  mu_x;
    theGeometry->muY         =  mu_y;
    theGeometry->N            =  N;
    theGeometry->joukowsky_x  = (double *)malloc(N * sizeof(double));
    theGeometry->joukowsky_y  = (double *)malloc(N * sizeof(double));
    theGeometry->dCircle1     = theGeometry->h * 2;
    theGeometry->dCircle2     = theGeometry->h * 2;
    theGeometry->dCircle3     = theGeometry->h * 2;
    theGeometry->hCircle1     = theGeometry->h * 0.5;
    theGeometry->hCircle2     = theGeometry->h * 0.5;
    theGeometry->hCircle3     = theGeometry->h * 0.25;
    
    geoMeshGenerate();
    geoMeshImport();
    geoSetDomainName(0,"Circle1");
    geoSetDomainName(1,"Foil");
    geoSetDomainName(2,"Circle2");
    geoSetDomainName(3,"Circle3");
    geoMeshWrite("../data/mesh.txt");

    printf("Mesh generated in %f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);
}
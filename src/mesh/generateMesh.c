#include "generateMesh.h"

void generateMesh(void){
    clock_t start = clock();
    double R = 1.2;
    double mu_x = -0.1;
    double mu_y = 0.1;
    int N = 100;
    

    geoInitialize();
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
    theGeometry->hCircle1     = theGeometry->h;
    theGeometry->hCircle2     = theGeometry->h;
    theGeometry->hCircle3     = theGeometry->h;

    theGeometry->elementType  = FEM_TRIANGLE;
    geoMeshGenerate();
    geoMeshImport();
    geoSetDomainName(0,"Circle1");
    geoSetDomainName(1,"Foil");
    geoSetDomainName(2,"Circle2");
    geoSetDomainName(3,"Circle3");
    geoMeshWrite("../../data/mesh.txt");
    geoFinalize();

    printf("Mesh generated in %f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);
}


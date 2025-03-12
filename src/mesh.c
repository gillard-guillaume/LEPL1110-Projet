#include "mesh.h"


double geoSize(double x, double y){
    return 0.1;
}

int generateSurface(double *x, double *y, int N) {
    int ierr;
    
    int *points = malloc((N+1) * sizeof(int));
    if (!points) {
        printf("Error: Could not allocate memory for points.\n");
        return 1;
    }
    
    for (int i = 0; i < N; i++) {
        points[i] = gmshModelOccAddPoint(x[i], y[i], 0, 0.1, -1, &ierr);
        if (ierr) {
            printf("Error: Could not add point %d.\n", i);
            return 1;
        }
    }
    points[N] = points[0];
    int splineID = gmshModelOccAddSpline(points, N+1, -1, NULL, 0, &ierr);
    int loopID = gmshModelOccAddCurveLoop(&splineID, 1, -1, &ierr);
    int surfaceID = gmshModelOccAddPlaneSurface(&loopID, 1, -1, &ierr);

    return surfaceID;
}

   

int getPoints(double *x, double *y, int n, char *filename){
    FILE *f = fopen(filename, "r");
    for (int i = 0; i < n; i++){
        fscanf(f, "%lf %lf\n", &x[i], &y[i]);
    }
    fclose(f);
    return 0;
}



int wing(){
    

    int ierr;

    gmshInitialize(0, NULL, 1, 0, &ierr);                         
    gmshModelAdd("Wing", &ierr);
    
    int N = 100;

    // Getting the points from the data files
    double *joukowsky_x = (double *)malloc(N * sizeof(double));
    double *joukowsky_y = (double *)malloc(N * sizeof(double));
    double *circle1_x = (double *)malloc(N * sizeof(double));
    double *circle1_y = (double *)malloc(N * sizeof(double));
    double *circle2_x = (double *)malloc(N * sizeof(double));
    double *circle2_y = (double *)malloc(N * sizeof(double));
    double *circle3_x = (double *)malloc(N * sizeof(double));
    double *circle3_y = (double *)malloc(N * sizeof(double));

    getPoints(joukowsky_x, joukowsky_y, N, "../data/joukowsky.dat");
    getPoints(circle1_x, circle1_y, N, "../data/circle1.dat");
    getPoints(circle2_x, circle2_y, N, "../data/circle2.dat");
    getPoints(circle3_x, circle3_y, N, "../data/circle3.dat");


    printf("🚀 Génération des surfaces...\n");
    int joukowsky = generateSurface(joukowsky_x, joukowsky_y, N);
    printf("✅ Surface de Joukowsky générée\n");
    int circle1 = generateSurface(circle1_x, circle1_y, N);
    printf("✅ Surface du cercle 1 générée\n");
    int circle2 = generateSurface(circle2_x, circle2_y, N);
    printf("✅ Surface du cercle 2 générée\n");
    int circle3 = generateSurface(circle3_x, circle3_y, N);
    printf("✅ Surface du cercle 3 générée\n");

    int joukowskyID[] = {2, joukowsky};
    int circle1ID[] = {2, circle1};
    int circle2ID[] = {2, circle2};
    int circle3ID[] = {2, circle3};

    gmshModelOccCut(joukowskyID,2,circle1ID ,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(joukowskyID,2,circle2ID ,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(joukowskyID,2,circle3ID ,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    gmshModelOccSynchronize(&ierr);
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);

    return 0;
}
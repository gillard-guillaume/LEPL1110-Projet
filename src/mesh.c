#include "mesh.h"

double interpolationH(double x, double x0, double x1, double h0, double h) {
    double q0 = pow((x-x1)/(x0-x1),2);
    double q1 = pow((x-x0)/(x1-x0),2);
    double P0 = h0 + (x-x0)*(-2/(x0-x1))*h0;
    double P1 = h + (x-x1)*(-2/(x1-x0))*h;
    return q0*P0 + q1*P1;
}


double geoSize(double x, double y){
    femGeo *theGeometry = geoGetGeometry();
    double h = theGeometry->h;

    double x0 = theGeometry->xCircle1;
    double y0 = theGeometry->yCircle1;
    double r0 = theGeometry->rCircle1;
    double d0 = theGeometry->dCircle1;
    double h0 = theGeometry->hCircle1;

    double x1 = theGeometry->xCircle2;
    double y1 = theGeometry->yCircle2;
    double r1 = theGeometry->rCircle2;
    double d1 = theGeometry->dCircle2;
    double h1 = theGeometry->hCircle3;

    double x2 = theGeometry->xCircle3;
    double y2 = theGeometry->yCircle3;
    double r2 = theGeometry->rCircle3;
    double d2 = theGeometry->dCircle3;
    double h2 = theGeometry->hCircle3;

    double hfinal = h;
    double d = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) - r0;
    if (d < d0) {
        double a = (-2*h + 2*h0)/(d0*d0*d0);
        double b = (3*h  - 3*h0)/(d0*d0);
        double c = 0;
        hfinal = a*d*d*d + b*d*d + c*d + h0; }
        
    d = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1)) - r1;
    if (d < d1) {
        double a = (-2*h + 2*h1)/(d1*d1*d1);
        double b = (3*h  - 3*h1)/(d1*d1);
        double c = 0;
        hfinal = fmin(hfinal,a*d*d*d + b*d*d + c*d + h1); }
    
    d = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2)) - r2;
    if (d < d2) {
        double a = (-2*h + 2*h2)/(d2*d2*d2);
        double b = (3*h  - 3*h2)/(d2*d2);
        double c = 0;
        hfinal = fmin(hfinal,a*d*d*d + b*d*d + c*d + h2); }

    return hfinal;
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

int wing(femGeo *theGeometry) {
    
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

    printf("Generating surfaces...\n");
    int joukowsky = generateSurface(joukowsky_x, joukowsky_y, N);
    printf("Joukowsky surface generated\n");
    int circle1 = generateSurface(circle1_x, circle1_y, N);
    printf("Circle 1 surface generated\n");
    int circle2 = generateSurface(circle2_x, circle2_y, N);
    printf("Circle 2 surface generated\n");
    int circle3 = generateSurface(circle3_x, circle3_y, N);
    printf("Circle 3 surface generated\n");

    int joukowskyID[] = {2, joukowsky};
    int circle1ID[] = {2, circle1};
    int circle2ID[] = {2, circle2};
    int circle3ID[] = {2, circle3};

    gmshModelOccCut(joukowskyID,2,circle1ID ,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(joukowskyID,2,circle2ID ,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(joukowskyID,2,circle3ID ,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    geoSetSizeCallback(geoSize);

    gmshModelOccSynchronize(&ierr);
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);

    return 0;
}
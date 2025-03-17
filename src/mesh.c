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

    // Joukowsky parameters
    double R = theGeometry->R;
    double muX = theGeometry->muX;
    double muY = theGeometry->muY;
    double hJoukowsky = theGeometry->hJoukowsky;
    
    // Circles parameters
    double xCircle1 = theGeometry->xCircle1, yCircle1 = theGeometry->yCircle1, rCircle1 = theGeometry->rCircle1, hCircle1 = theGeometry->hCircle1, dCircle1 = theGeometry->dCircle1;
    double xCircle2 = theGeometry->xCircle2, yCircle2 = theGeometry->yCircle2, rCircle2 = theGeometry->rCircle2, hCircle2 = theGeometry->hCircle2, dCircle2 = theGeometry->dCircle2;
    double xCircle3 = theGeometry->xCircle3, yCircle3 = theGeometry->yCircle3, rCircle3 = theGeometry->rCircle3, hCircle3 = theGeometry->hCircle3, dCircle3 = theGeometry->dCircle3;

    // Getting back to the Joukowsky plane
    double complex zeta = inverseJoukowskyTransform(x + I*y);
    double xJoukowsky = creal(zeta);
    double yJoukowsky = cimag(zeta);
    // Computing the distance to the center of the joukowski circle
    double dJoukowsky = sqrt(pow(xJoukowsky - muX, 2) + pow(yJoukowsky - muY, 2));
    // Computing the distance to the circles
    dCircle1 = sqrt(pow(x - xCircle1, 2) + pow(y - yCircle1, 2));
    dCircle2 = sqrt(pow(x - xCircle2, 2) + pow(y - yCircle2, 2));
    dCircle3 = sqrt(pow(x - xCircle3, 2) + pow(y - yCircle3, 2));

    double hres = h;
    double hres0 = h;
    double hres1 = h, hres2 = h, hres3 = h;

    if (dJoukowsky < 0) return hres;
    if (dCircle1 < 0) return hres1;
    if (dCircle2 < 0) return hres2;
    if (dCircle3 < 0) return hres3;

    if (dJoukowsky > 0.75*R) return interpolationH(dJoukowsky, R, 0.9*R, hJoukowsky, h);
    if (dCircle1 < 1.25*rCircle1) return interpolationH(dCircle1, rCircle1, 1.1*rCircle1, hCircle1, h);
    if (dCircle2 < 1.25*rCircle2) return interpolationH(dCircle2, rCircle2, 1.1*rCircle2, hCircle2, h);
    if (dCircle3 < 1.25*rCircle3) return interpolationH(dCircle3, rCircle3, 1.1*rCircle3, hCircle3, h);


    return hres;

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
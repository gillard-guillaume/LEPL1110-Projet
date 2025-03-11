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
    int splineID = gmshModelOccAddSpline(points, N, -1, NULL, 0, &ierr);
    if (ierr) {
        printf("Error: Could not add spline.\n");
        free(points);
        return 1;
    }

    free(points);

    int loopID = gmshModelOccAddCurveLoop(&splineID, 1, -1, &ierr);
    if (ierr) {
        printf("Error: Could not add curve loop.\n");
        return 1;
    }   

    int surfaceID = gmshModelOccAddPlaneSurface(&loopID, 1, -1, &ierr);
    if (ierr) {
        printf("Error: Could not add plane surface.\n");
        return 1;
    }

    return surfaceID;
}


int getPoints (const char *filename, double **x, double **y, int *N){
    char path[256];
    sprintf(path, "../data/%s", filename);
    FILE *file = fopen(path, "r");
    if (!file) {
        printf("Error: Unable to open file %s\n", filename);
        return -1;
    }

    
    int count = 0;
    double temp_x, temp_y;
    

    while (fscanf(file, "%lf %lf", &temp_x, &temp_y) == 2) {
        count++;
    }
    rewind(file);  

    *x = (double *)malloc(count * sizeof(double));
    *y = (double *)malloc(count * sizeof(double));
    if (!(*x) || !(*y)) {
        printf("Error: Memory allocation failed.\n");
        fclose(file);
        return -1;
    }

    for (int i = 0; i < count; i++) {
        fscanf(file, "%lf %lf", &((*x)[i]), &((*y)[i]));
    }

    fclose(file);
    *N = count;  
    return 0;
}



int wing(){
    

    int ierr;

    gmshInitialize(0, NULL, 1, 0, &ierr);                         
    gmshModelAdd("Wing", &ierr);
    
    // Getting the points from the data files
    double *joukowsky_x, *joukowsky_y, *circle1_x, *circle1_y, *circle2_x, *circle2_y, *circle3_x, *circle3_y;
    int N;
    if (getPoints("joukowsky.dat", &joukowsky_x, &joukowsky_y, &N) || getPoints("circle1.dat", &circle1_x, &circle1_y, &N) || getPoints("circle2.dat", &circle2_x, &circle2_y, &N) || getPoints("circle3.dat", &circle3_x, &circle3_y, &N)) {
        return 1;
    }
    printf("🚀 Génération des surfaces...\n");
    int joukowsky = generateSurface(joukowsky_x, joukowsky_y, N);
    int circle1 = generateSurface(circle1_x, circle1_y, N);
    int circle2 = generateSurface(circle2_x, circle2_y, N);
    int circle3 = generateSurface(circle3_x, circle3_y, N);

    int objectDimTags[] = {2, joukowsky};
    int toolDimTags[] = {2, circle1, 2, circle2, 2, circle3};

    int *outDimTags = NULL;
    size_t outDimTags_n = 0;

    gmshModelOccCut(objectDimTags, 2, toolDimTags, 6, &outDimTags, &outDimTags_n, NULL, NULL, NULL, -1, 1, 1, &ierr);
    if (ierr) {
        printf("Error: Could not cut surfaces.\n");
        return 1;
    }

    int wing = outDimTags[1];

    gmshModelOccSynchronize(&ierr);
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);

    return 0;
}
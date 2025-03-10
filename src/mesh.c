#include "mesh.h"


int wing_surface(double *x, double *y, int N) {
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
    if (ierr) {
        printf("Error: Could not add spline.\n");
        return 1;
    }

    int loopID = gmshModelOccAddCurveLoop([splineID], 1, &ierr);
    if (ierr) {
        printf("Error: Could not add curve loop.\n");
        return 1;
    }   

    int surfaceID = gmshModelOccAddPlaneSurface([loopID], 1, &ierr);
    if (ierr) {
        printf("Error: Could not add plane surface.\n");
        return 1;
    }

    int wing[] = {2, surfaceID};
}
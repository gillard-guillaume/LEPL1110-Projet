#include "joukowsky.h"


double complex joukowskyTransform(double complex zeta){
    // Joukowsky transformation
    //z = zeta + 1/zeta;
    double complex xi = creal(zeta);
    double complex eta = cimag(zeta);
    double x = xi * (1 + 1/(pow(xi, 2) + pow(eta, 2)));
    double y = eta * (1 - 1/(pow(xi, 2) + pow(eta, 2)));
    return CMPLX(x, y);
}

double complex inverseJoukowskyTransform(double complex z){
    // Inverse Joukowsky transformation
    //zeta = 0.5*(z + sqrt(z^2 - 4))
    double xi = creal(z);
    double eta = cimag(z);
    double complex zeta = 0.5 * (z + csqrt(pow(z, 2) - 4));
    return zeta;
}


int joukowsky(femGeo *theGeometry){

    FILE *file = fopen("../data/joukowsky.dat", "w");
    if (!file) {
        printf("Error: Could not open file for writing0.\n");
        return 1;
    }

    double R = theGeometry->R;
    double mu_x = theGeometry->muX;
    double mu_y = theGeometry->muY;
    double N = theGeometry->N;
    double *joukowsky_x = theGeometry->joukowsky_x;
    double *joukowsky_y = theGeometry->joukowsky_y;

    // Generating the airfoil points
    double x_min = 0.0;
    double x_max = 0.0;
    for (int i = 0; i < N; i++){
        double theta = 2.0 * M_PI * i / N;
        double complex zeta = R * cos(theta) + I * R * sin(theta) + mu_x + I * mu_y;
        double complex z = joukowskyTransform(zeta);
        double x = creal(z);
        double y = cimag(z);
        
        // Finding the min and max x values
        if (x < x_min) x_min = x;
        if (x > x_max) x_max = x;
        joukowsky_x[i] = x;
        joukowsky_y[i] = y;
        // Writing the points to the file
        fprintf(file, "%f %f\n", creal(z), cimag(z));
    }
    fclose(file);

    // Finding insides circles centers and radius
    double wing_length = fabs(x_max) + fabs(x_min);
    theGeometry->Length = wing_length;

    double x_center1 = 0.0 - wing_length/3;
    double x_center2 = 0.0;
    double x_center3 = 0.0 + wing_length/3;

    theGeometry->xCircle1 = x_center1;
    theGeometry->xCircle2 = x_center2;
    theGeometry->xCircle3 = x_center3;

    double x, y;
    double y_max1 = 0.0, y_min1 = 0.0;
    double y_max2 = 0.0, y_min2 = 0.0;
    double y_max3 = 0.0, y_min3 = 0.0;


    // Reading the file to find y_max and y_min at x_center1, x_center2, x_center3
    for (int i = 0; i < N; i++){
        x = joukowsky_x[i];
        y = joukowsky_y[i];
        if (fabs(x - x_center1) < 0.05) {  
            if (y > y_max1) y_max1 = y;
            if (y < y_min1) y_min1 = y;
        }
        if (fabs(x - x_center2) < 0.05) {
            if (y > y_max2) y_max2 = y;
            if (y < y_min2) y_min2 = y;
        }
        if (fabs(x - x_center3) < 0.05) { 
            if (y > y_max3) y_max3 = y;
            if (y < y_min3) y_min3 = y;
        }
    }

    double y_center1 = (y_max1 + y_min1) / 2;
    double y_center2 = (y_max2 + y_min2) / 2;
    double y_center3 = (y_max3 + y_min3) / 2;

    theGeometry->yCircle1 = y_center1;
    theGeometry->yCircle2 = y_center2;
    theGeometry->yCircle3 = y_center3;

    
    double height1 = fabs(y_max1) + fabs(y_min1);
    double height2 = fabs(y_max2) + fabs(y_min2);
    double height3 = fabs(y_max3) + fabs(y_min3);
    double c1 = 0.3;
    double c2 = 0.3;
    double c3 = 0.3;
    double R1 = c1 * height1;
    double R2 = c2 * height2;
    double R3 = c3 * height3;

    theGeometry->rCircle1 = R1;
    theGeometry->rCircle2 = R2;
    theGeometry->rCircle3 = R3;

    return 0;
}
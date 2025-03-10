// Joukowsky c implementation

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif




double complex joukowski(double complex zeta){
    // Joukowsky transformation
    //z = zeta + 1/zeta;
    double complex xi = creal(zeta);
    double complex eta = cimag(zeta);
    double x = xi * (1 + 1/(pow(xi, 2) + pow(eta, 2)));
    double y = eta * (1 - 1/(pow(xi, 2) + pow(eta, 2)));
    return CMPLX(x, y);
}


int main(int argc, char *argv[]){
    if (argc != 5){
        printf("Usage: %s <R> <mu_y> <mu_x> <N>\n", argv[0]);
        printf("  R     = Circle radius (should be >1)\n");
        printf("  mu_x  = X offset of the circle center\n");
        printf("  mu_y  = Y offset of the circle center\n");
        printf("  N     = Number of points on the circle\n");
        return 1;
    }

    // Reading input arguments
    double R = atof(argv[1]);
    double mu_x = atof(argv[2]);
    double mu_y = atof(argv[3]);
    int N = atoi(argv[4]);

    // Checking the input arguments
    if (R <= 1.0 || N < 3) {
        printf("Error: R must be >1 and N must be at least 3.\n"); // With R=1 for some angle the transformation will fail
        return 1;
    }


    FILE *file = fopen("../data/joukowsky.dat", "w");
    if (!file) {
        printf("Error: Could not open file for writing0.\n");
        return 1;
    }


    // Generating the airfoil points
    double x_min = 0.0;
    double x_max = 0.0;
    for (int i = 0; i <= N; i++){
        double theta = 2.0 * M_PI * i / N;
        double complex zeta = R * cos(theta) + I * R * sin(theta) + mu_x + I * mu_y;
        double complex z = joukowski(zeta);
        double x = creal(z);
        double y = cimag(z);
        
        // Finding the min and max x values
        if (x < x_min) x_min = x;
        if (x > x_max) x_max = x;
        fprintf(file, "%f %f\n", creal(z), cimag(z));
    }

    fclose(file);

    // Generating insides circles 
    double wing_length = fabs(x_max) + fabs(x_min);

    double x_center1 = 0.0 - wing_length/3;
    double x_center2 = 0.0;
    double x_center3 = 0.0 + wing_length/3;


    file = fopen("../data/joukowsky.dat", "r");
    if (!file) {
        printf("Error: Could not open file for reading.\n");
        return 1;
    }
    double x, y;

    double y_max1 = 0.0, y_min1 = 0.0;
    double y_max2 = 0.0, y_min2 = 0.0;
    double y_max3 = 0.0, y_min3 = 0.0;


    // Reading the file to find y_max and y_min at x_center1, x_center2, x_center3
    while (fscanf(file, "%lf %lf", &x, &y) == 2) {
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

    fclose(file);

    double y_center1 = (y_max1 + y_min1) / 2;
    double y_center2 = (y_max2 + y_min2) / 2;
    double y_center3 = (y_max3 + y_min3) / 2;

    
    double height1 = fabs(y_max1) + fabs(y_min1);
    double height2 = fabs(y_max2) + fabs(y_min2);
    double height3 = fabs(y_max3) + fabs(y_min3);
    double c1 = 0.3;
    double c2 = 0.3;
    double c3 = 0.3;
    double R1 = c1 * height1;
    double R2 = c2 * height2;
    double R3 = c3 * height3;


    file = fopen("../data/circle1.dat", "w");
    if (!file) {
        printf("Error: Could not open file for writing1.\n");
        return 1;
    }

    for (int i = 0; i <= N; i++){
        double theta = 2.0 * M_PI * i / N;
        double complex z1 = x_center1 + R1 * cos(theta) + I * (y_center1 + R1 * sin(theta));
        fprintf(file, "%f %f\n", creal(z1), cimag(z1));
    }
    fclose(file);


    file = fopen("../data/circle2.dat", "w");
    if (!file) {
        printf("Error: Could not open file for writing.\n");
        return 1;
    }

    for (int i = 0; i <= N; i++){
        double theta = 2.0 * M_PI * i / N;
        double complex z2 = x_center2 + R2 * cos(theta) + I * (y_center2 + R2 * sin(theta));
        fprintf(file, "%f %f\n", creal(z2), cimag(z2));
    }
    fclose(file);

    file = fopen("../data/circle3.dat", "w");
    if (!file) {
        printf("Error: Could not open file for writing.\n");
        return 1;
    }
    
    for (int i = 0; i <= N; i++){
        double theta = 2.0 * M_PI * i / N;
        double complex z3 = x_center3 + R3 * cos(theta) + I * (y_center3 + R3 * sin(theta));
        fprintf(file, "%f %f\n", creal(z3), cimag(z3));
    }
 
    fclose(file);


}
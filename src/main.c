#include <main.h>

int main(int argc, char *argv[]) {
    int ierr;
    
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
        printf("Error: R must be >1 and N must be at least 3.\n"); // With R=1 for some angle the joukowsky transformation will fail
        return 1;
    }

    // Generating Joukowsky airfoil and circles points
    printf("Generating Joukowsky airfoil and circles points...\n");
    if (joukowsky(R, mu_x, mu_y, N) != 0) {
        printf("Error: Joukowsky airfoil and circles generation failed.\n");
        return 1;
    }
    printf("Joukowsky airfoil and circles points generated successfully.\n");


    printf("Generating wing mesh...\n");
    if (wing() != 0) {
        printf("Error: Wing generation failed.\n");
        return 1;
    }
    printf("Wing mesh generated successfully.\n");

    // Save the mesh to a file
    gmshWrite("../data/wing.msh", &ierr);
    if (ierr) {
        printf("Error: Could not write mesh file.\n");
        return 1;
    }

    printf("Mesh saved to '../data/wing.msh'\n");

    // Display the mesh using Gmsh GUI
    gmshFltkRun(&ierr);

    // Finalize Gmsh
    gmshFinalize(&ierr);
    
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "../gmsh/gmsh-4.13.1-Linux64-sdk/include/gmshc.h"

int main() {
    int ierr;
    
    printf("🚀 Generating wing mesh...\n");
    if (wing() != 0) {
        printf("Error: Wing generation failed.\n");
        return 1;
    }
    printf("✅ Wing mesh generated successfully.\n");

    // Save the mesh to a file
    gmshWrite("../data/wing.msh", &ierr);
    if (ierr) {
        printf("Error: Could not write mesh file.\n");
        return 1;
    }

    printf("✅ Mesh saved to '../data/wing.msh'\n");

    // Display the mesh using Gmsh GUI
    gmshFltkRun(&ierr);

    // Finalize Gmsh
    gmshFinalize(&ierr);
    
    return 0;
}

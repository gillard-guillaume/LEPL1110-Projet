#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "../gmsh/gmsh-4.13.1-Linux64-sdk/include/gmshc.h"

int main() {
    int ierr;
    
    printf("🚀 Génération du maillage de l'aile...\n");
    if (wing() != 0) {
        printf("Error: Wing generation failed.\n");
        return 1;
    }
    printf("✅ Maillage de l'aile généré\n");

    gmshWrite("../data/wing.msh", &ierr);
    if (ierr) {
        printf("Error: Could not write mesh file.\n");
        return 1;
    }

    printf("✅ Maillage généré et sauvegardé dans '../data/wing.msh'\n");

    gmshFinalize(&ierr);
    
    return 0;
}

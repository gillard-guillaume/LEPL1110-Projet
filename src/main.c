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

    geoInitialize();
    femGeo *theGeometry = geoGetGeometry();
    theGeometry->R = R;
    theGeometry->muX = mu_x;
    theGeometry->muY = mu_y;
    theGeometry->h = 0.1;
    theGeometry->elementType = FEM_TRIANGLE;
    theGeometry->geoSize = geoSize;


    theGeometry->h = 2.5;
    theGeometry->hCircle1 = 0.6;
    theGeometry->hCircle2 = 2.2;
    theGeometry->hCircle3 = 0.3;

    // Generating Joukowsky airfoil and circles points
    printf("Generating Joukowsky airfoil and circles points...\n");
    if (joukowsky(R, mu_x, mu_y, N, theGeometry) != 0) {
        printf("Error: Joukowsky airfoil and circles generation failed.\n");
        return 1;
    }
    printf("Joukowsky airfoil and circles points generated successfully.\n");


    printf("Generating wing mesh...\n");
    if (wing(theGeometry) != 0) {
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


    geoMeshImport();

    gmshFltkInitialize(&ierr);
    gmshFltkRun(&ierr);

//
//  -3- Champ de la taille de r�f�rence du maillage
//

    double *meshSizeField = malloc(theGeometry->theNodes->nNodes*sizeof(double));
    femNodes *theNodes = theGeometry->theNodes;
    for(int i=0; i < theNodes->nNodes; ++i)
        meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i]);
    double hMin = femMin(meshSizeField,theNodes->nNodes);  
    double hMax = femMax(meshSizeField,theNodes->nNodes);  
    printf(" ==== Global requested h : %14.7e \n",theGeometry->h);
    printf(" ==== Minimum h          : %14.7e \n",hMin);
    printf(" ==== Maximum h          : %14.7e \n",hMax);
 
//
//  -4- Visualisation du maillage
//  
    
    int mode = 1; // Change mode by pressing "j", "k", "l"
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[256];
    double pos[2] = {20,460};
 
 
    GLFWwindow* window = glfemInit("EPL1110 : Mesh generation ");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
    
    
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
    //    glfemChangeState(&mode, theMeshes->nMesh);
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}

        
        if (t-told > 0.5) {freezingButton = FALSE; }
            
        
        
         
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements, meshSizeField);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
 
            
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
 
            
            
            }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            
            
            
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
 
            
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
            }
            
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(meshSizeField);  
    geoFinalize();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}

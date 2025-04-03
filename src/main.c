/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Elasticite lineaire plane
 *
 *  Copyright (C) 2024 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "main.h"

int main(void)
{  
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n\n\n");

    femGeo* theGeometry = geoGetGeometry();
    generateMesh();
    geoMeshRead("../data/mesh.txt");

            
//
//  -2- Creation probleme 
//
    // Aluminium 7075-T6  
    double E   = 71.7e9;
    double nu  = 0.33;
    double rho = 2.810e3; 
    double g   = 9.81 * 1;
    femRenumType renumType = FEM_XNUM;
    femSolverType solverType = FEM_CHOV;
    theGeometry->elementType  = FEM_TRIANGLE;
    femElasticCase caseType = PLANAR_STRAIN;


    clock_t start = clock();
    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g, caseType, renumType);
    femElasticityAddBoundaryCondition(theProblem, "Foil", NEUMANN_Y, 0);
    femElasticityAddBoundaryCondition(theProblem,"Circle1",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Circle1",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Circle2",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Circle2",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Circle3",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Circle3",DIRICHLET_Y,0.0);

    femElasticityPrint(theProblem);
    printf("Problem created in %f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);
    start = clock();
    double *theSoluce = femElasticitySolve(theProblem, solverType);
    printf("Problem solved in %f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);
   
//
//  -3- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//
    
    femNodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 1e5;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    
    for (int i=0; i<theNodes->nNodes; i++){
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                   theSoluce[2*i+1]*theSoluce[2*i+1]); }
  
    double hMin = femMin(normDisplacement,theNodes->nNodes);  
    double hMax = femMax(normDisplacement,theNodes->nNodes);  
    printf(" ==== Minimum displacement          : %14.7e \n",hMin);
    printf(" ==== Maximum displacement          : %14.7e \n",hMax);
 
//
//  -4- Visualisation du maillage
//  
    
    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];
   
 
    GLFWwindow* window = glfemInit("EPL1110 : Linear elasticity ");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'S') == GLFW_PRESS) { mode = 2;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        
        if (t-told > 0.5) {freezingButton = FALSE; }
        if (mode == 2){
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->system,theProblem->system->size,w,h);
        }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
             glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);  }
            
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(normDisplacement);
    femElasticityFree(theProblem) ; 
    geoFinalize();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}

 

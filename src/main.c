/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Elasticite lineaire plane
 *
 *  Copyright (C) 2024 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include <main.h>



int main(void)
{  
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n\n\n");

    double R = 1.2;
    double mu_x = -0.1;
    double mu_y = 0.1;
    int N = 100;
       
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();    
    theGeometry->h            =  0.1;  
    theGeometry->R            =  R;
    theGeometry->muX         =  mu_x;
    theGeometry->muY         =  mu_y;
    theGeometry->N            =  N;
    theGeometry->joukowsky_x  = (double *)malloc(N * sizeof(double));
    theGeometry->joukowsky_y  = (double *)malloc(N * sizeof(double));
    theGeometry->dCircle1     = theGeometry->h * 2;
    theGeometry->dCircle2     = theGeometry->h * 2;
    theGeometry->dCircle3     = theGeometry->h * 2;
    theGeometry->hCircle1     = theGeometry->h;
    theGeometry->hCircle2     = theGeometry->h;
    theGeometry->hCircle3     = theGeometry->h;

    theGeometry->elementType  = FEM_TRIANGLE;
   
    clock_t start = clock();

// geoMeshGenerate();
// geoMeshImport();
// geoSetDomainName(0,"Circle1");
// geoSetDomainName(1,"Foil");
// geoSetDomainName(2,"Circle2");
    // geoSetDomainName(3,"Circle3");
    // geoMeshWrite("../data/mesh.txt");
    geoMeshRead("../data/mesh.txt");

    printf("Mesh generated in %f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);
            
//
//  -2- Creation probleme 
//
    // Aluminium 7075-T6  
    double E   = 71.7e9;
    double nu  = 0.33;
    double rho = 2.810e3; 
    double g   = 9.81 * 1;

    start = clock();
    femRenumType renumType = FEM_XNUM;
    femSolverType solverType = FEM_CHOV;

    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN, renumType);
    femElasticityAddBoundaryCondition(theProblem, "Foil", NEUMANN_Y, 1e4);
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
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        
        if (t-told > 0.5) {freezingButton = FALSE; }
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
 
  

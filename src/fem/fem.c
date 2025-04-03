/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

 #include "fem.h"


 static const double _gaussQuad4Xsi[4]    = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626};
 static const double _gaussQuad4Eta[4]    = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626};
 static const double _gaussQuad4Weight[4] = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
 static const double _gaussTri3Xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
 static const double _gaussTri3Eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
 static const double _gaussTri3Weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};
 static const double _gaussEdge2Xsi[2]    = { 0.577350269189626,-0.577350269189626};
 static const double _gaussEdge2Weight[2] = { 1.000000000000000, 1.000000000000000};
 double **A_copy = NULL;
 double *B_copy  = NULL;
 
 femIntegration *femIntegrationCreate(int n, femElementType type)
 {
     femIntegration *theRule = malloc(sizeof(femIntegration));
     if (type == FEM_QUAD && n == 4) {
         theRule->n      = 4;
         theRule->xsi    = _gaussQuad4Xsi;
         theRule->eta    = _gaussQuad4Eta;
         theRule->weight = _gaussQuad4Weight; }
     else if (type == FEM_TRIANGLE && n == 3) {
         theRule->n      = 3;
         theRule->xsi    = _gaussTri3Xsi;
         theRule->eta    = _gaussTri3Eta;
         theRule->weight = _gaussTri3Weight; }
     else if (type == FEM_EDGE && n == 2) {
         theRule->n      = 2;
         theRule->xsi    = _gaussEdge2Xsi;
         theRule->eta    = NULL;
         theRule->weight = _gaussEdge2Weight; }
     else Error("Cannot create such an integration rule !");
     return theRule; 
 }
 
 void femIntegrationFree(femIntegration *theRule)
 {
     free(theRule);
 }
 
 void _q1c0_x(double *xsi, double *eta) 
 {
     xsi[0] =  1.0;  eta[0] =  1.0;
     xsi[1] = -1.0;  eta[1] =  1.0;
     xsi[2] = -1.0;  eta[2] = -1.0;
     xsi[3] =  1.0;  eta[3] = -1.0;
 }
 
 void _q1c0_phi(double xsi, double eta, double *phi)
 {
     phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;  
     phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
     phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
     phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
 }
 
 void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
 {
     dphidxsi[0] =   (1.0 + eta) / 4.0;  
     dphidxsi[1] = - (1.0 + eta) / 4.0;
     dphidxsi[2] = - (1.0 - eta) / 4.0;
     dphidxsi[3] =   (1.0 - eta) / 4.0;
     dphideta[0] =   (1.0 + xsi) / 4.0;  
     dphideta[1] =   (1.0 - xsi) / 4.0;
     dphideta[2] = - (1.0 - xsi) / 4.0;
     dphideta[3] = - (1.0 + xsi) / 4.0;
 
 }
 
 void _p1c0_x(double *xsi, double *eta) 
 {
     xsi[0] =  0.0;  eta[0] =  0.0;
     xsi[1] =  1.0;  eta[1] =  0.0;
     xsi[2] =  0.0;  eta[2] =  1.0;
 }
 
 void _p1c0_phi(double xsi, double eta, double *phi)
 {
     phi[0] = 1 - xsi - eta;  
     phi[1] = xsi;
     phi[2] = eta;
 }
 
 void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
 {
     dphidxsi[0] = -1.0;  
     dphidxsi[1] =  1.0;
     dphidxsi[2] =  0.0;
     dphideta[0] = -1.0;  
     dphideta[1] =  0.0;
     dphideta[2] =  1.0;
 }
 void _e1c0_x(double *xsi) 
 {
     xsi[0] = -1.0;  
     xsi[1] =  1.0;  
 }
 
 void _e1c0_phi(double xsi,  double *phi)
 {
     phi[0] = (1 - xsi) / 2.0;  
     phi[1] = (1 + xsi) / 2.0;
 }
 
 void _e1c0_dphidx(double xsi, double *dphidxsi)
 {
     dphidxsi[0] = -0.5;  
     dphidxsi[1] =  0.5;
 }
 
 
 femDiscrete *femDiscreteCreate(int n, femElementType type)
 {
     femDiscrete *theSpace = malloc(sizeof(femDiscrete));
     if (type == FEM_QUAD && n == 4) {
         theSpace->n       = 4;
         theSpace->x2      = _q1c0_x;
         theSpace->phi2    = _q1c0_phi;
         theSpace->dphi2dx = _q1c0_dphidx; }
     else if (type == FEM_TRIANGLE && n == 3) {
         theSpace->n       = 3;
         theSpace->x2      = _p1c0_x;
         theSpace->phi2    = _p1c0_phi;
         theSpace->dphi2dx = _p1c0_dphidx; }
     else if (type == FEM_EDGE && n == 2) {
         theSpace->n       = 2;
         theSpace->x       = _e1c0_x;
         theSpace->phi     = _e1c0_phi;
         theSpace->dphidx  = _e1c0_dphidx; }
     else Error("Cannot create such a discrete space !");
     return theSpace; 
 }
 
 void femDiscreteFree(femDiscrete *theSpace)
 {
     free(theSpace);
 }
 
 void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
 {
     mySpace->x2(xsi,eta);
 }
 
 void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
 {
     mySpace->phi2(xsi,eta,phi);
 }
 
 void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
 {
     mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
 }
 
 void femDiscreteXsi(femDiscrete* mySpace, double *xsi)
 {
     mySpace->x(xsi);
 }
 
 void femDiscretePhi(femDiscrete* mySpace, double xsi, double *phi)
 {
     mySpace->phi(xsi,phi);
 }
 
 void femDiscreteDphi(femDiscrete* mySpace, double xsi, double *dphidxsi)
 {
     mySpace->dphidx(xsi,dphidxsi);
 }
 
 void femDiscretePrint(femDiscrete *mySpace)
 {
     int i,j;
     int n = mySpace->n;
     double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];
     
     femDiscreteXsi2(mySpace,xsi,eta);
     for (i=0; i < n; i++) {
         
         femDiscretePhi2(mySpace,xsi[i],eta[i],phi);
         femDiscreteDphi2(mySpace,xsi[i],eta[i],dphidxsi,dphideta);
 
         for (j=0; j < n; j++)  {
             printf("(xsi=%+.1f,eta=%+.1f) : ",xsi[i],eta[i]);
             printf(" phi(%d)=%+.1f",j,phi[j]);  
             printf("   dphidxsi(%d)=%+.1f",j,dphidxsi[j]);  
             printf("   dphideta(%d)=%+.1f \n",j,dphideta[j]);  }
         printf(" \n"); }
 }
 
 
 femFullSystem *femFullSystemCreate(int size)
 {
     femFullSystem *theSystem = malloc(sizeof(femFullSystem));
     femFullSystemAlloc(theSystem, size);
     femFullSystemInit(theSystem);
 
     return theSystem; 
 }
 
 void femFullSystemFree(femFullSystem *theSystem)
 {
     free(theSystem->A);
     free(theSystem->B);
     free(theSystem);
 }
 
 void femFullSystemAlloc(femFullSystem *mySystem, int size)
 {
     int i;  
     double *elem = malloc(sizeof(double) * size * (size+1)); 
     mySystem->A = malloc(sizeof(double*) * size); 
     mySystem->B = elem;
     mySystem->A[0] = elem + size;  
     mySystem->size = size;
     for (i=1 ; i < size ; i++) 
         mySystem->A[i] = mySystem->A[i-1] + size;
 }
 
 void femFullSystemInit(femFullSystem *mySystem)
 {
     int i,size = mySystem->size;
     for (i=0 ; i < size*(size+1) ; i++) 
         mySystem->B[i] = 0;}
 
 
 void femFullSystemPrint(femFullSystem *mySystem)
 {
     double  **A, *B;
     int     i, j, size;
     
     A    = mySystem->A;
     B    = mySystem->B;
     size = mySystem->size;
     
     for (i=0; i < size; i++) {
         for (j=0; j < size; j++)
             if (A[i][j] == 0)  printf("         ");   
             else               printf(" %+.1e",A[i][j]);
         printf(" :  %+.1e \n",B[i]); }
 }
 
 double* femFullSystemEliminate(femFullSystem *mySystem, femSolverType solver)
 {
     double  **A, *B, factor;
     int     i, j, k, size;
     
     A    = mySystem->A;
     B    = mySystem->B;
     size = mySystem->size;
 
     /* Gauss elimination */
 
     for (k=0; k < size; k++) {
         if ( fabs(A[k][k]) <= 1e-16 ) {
             printf("Pivot index %d  ",k);
             printf("Pivot value %e  ",A[k][k]);
             Error("Cannot eliminate with such a pivot"); }
 
         for (i = k+1 ; i <  size; i++) {
             factor = A[i][k] / A[k][k];
             for (j = k+1 ; j < size; j++) 
                 A[i][j] = A[i][j] - A[k][j] * factor;
             B[i] = B[i] - B[k] * factor; 
         }
     }
     
     /* Back-substitution */
     
     for (i = size-1; i >= 0 ; i--) {
         factor = 0;
         for (j = i+1 ; j < size; j++)
             factor += A[i][j] * B[j];
         B[i] = ( B[i] - factor)/A[i][i]; }
 
     return(mySystem->B);    
 }
 
 void  femFullSystemConstrain(femFullSystem *mySystem, 
                              int myNode, double myValue) 
 {
     double  **A, *B;
     int     i, size;
     
     A    = mySystem->A;
     B    = mySystem->B;
     size = mySystem->size;
     
     for (i=0; i < size; i++) {
         B[i] -= myValue * A[i][myNode];
         A[i][myNode] = 0; }
     
     for (i=0; i < size; i++) 
         A[myNode][i] = 0; 
     
     A[myNode][myNode] = 1;
     B[myNode] = myValue;
 }
 
 
 femProblem *femElasticityCreate(femGeo* theGeometry, 
                   double E, double nu, double rho, double g, 
                   femElasticCase iCase, femRenumType renumType)
 {
     femProblem *theProblem = malloc(sizeof(femProblem));
     theProblem->E   = E;
     theProblem->nu  = nu;
     theProblem->g   = g;
     theProblem->rho = rho;
     
     if (iCase == PLANAR_STRESS) {
         theProblem->A = E/(1-nu*nu);
         theProblem->B = E*nu/(1-nu*nu);
         theProblem->C = E/(2*(1+nu)); }
     else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
         theProblem->A = E*(1-nu)/((1+nu)*(1-2*nu));
         theProblem->B = E*nu/((1+nu)*(1-2*nu));
         theProblem->C = E/(2*(1+nu)); }
 
     theProblem->planarStrainStress = iCase;
     theProblem->nBoundaryConditions = 0;
     theProblem->conditions = NULL;
     
     int size = 2*theGeometry->theNodes->nNodes;
     theProblem->constrainedNodes = malloc(size*sizeof(int));
     for (int i=0; i < size; i++) 
         theProblem->constrainedNodes[i] = -1;
     theProblem->geometry = theGeometry;  
     if (theGeometry->theElements->nLocalNode == 3) {
         theProblem->space    = femDiscreteCreate(3,FEM_TRIANGLE);
         theProblem->rule     = femIntegrationCreate(3,FEM_TRIANGLE);}
     if (theGeometry->theElements->nLocalNode == 4) {
         theProblem->space    = femDiscreteCreate(4,FEM_QUAD);
         theProblem->rule     = femIntegrationCreate(4,FEM_QUAD); }
 
     theProblem->spaceEdge    = femDiscreteCreate(2,FEM_EDGE);
     theProblem->ruleEdge     = femIntegrationCreate(2,FEM_EDGE);
     theProblem->system   = femFullSystemCreate(size);
     
     femMesh *theMesh = theGeometry->theElements;
     for (int i=0; i < 10; i++) printf("%d ",theMesh->nodes->number[i]);
     printf("\n");
     // femMeshRenumber(theMesh,renumType);
     for (int i=0; i < 10; i++) printf("%d ",theMesh->nodes->number[i]);
 
     return theProblem;
 }
 
 void femElasticityFree(femProblem *theProblem)
 {
     femFullSystemFree(theProblem->system);
     femIntegrationFree(theProblem->rule);
     femDiscreteFree(theProblem->space);
     free(theProblem->geometry->theNodes->number);
     free(theProblem->conditions);
     free(theProblem->constrainedNodes);
     free(theProblem);
 }
     
 void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value,  double (*profile)(double x, double y))
 {
     int iDomain = geoGetDomain(nameDomain);
     if (iDomain == -1)  Error("Undefined domain :-(");
 
     femBoundaryCondition* theBoundary = malloc(sizeof(femBoundaryCondition));
     theBoundary->domain = theProblem->geometry->theDomains[iDomain];
     theBoundary->value = value;
     theBoundary->type = type;
     theBoundary->profile =  profile;
     theProblem->nBoundaryConditions++;
     int size = theProblem->nBoundaryConditions;
     
     if (theProblem->conditions == NULL)
         theProblem->conditions = malloc(size*sizeof(femBoundaryCondition*));
     else 
         theProblem->conditions = realloc(theProblem->conditions, size*sizeof(femBoundaryCondition*));
     theProblem->conditions[size-1] = theBoundary;
     
     int shift=-1;
     if (type == DIRICHLET_X)  shift = 0;      
     if (type == DIRICHLET_Y)  shift = 1;  
     if (shift == -1) return; 
     int *elem = theBoundary->domain->elem;
     int nElem = theBoundary->domain->nElem;
     for (int e=0; e<nElem; e++) {
         for (int i=0; i<2; i++) {
             int node = theBoundary->domain->mesh->elem[2*elem[e]+i];
             theProblem->constrainedNodes[2*node+shift] = size-1; }}    
 }
 
 void femElasticityAssembleElements(femProblem *theProblem){
     femFullSystem  *theSystem   = theProblem->system;
     femIntegration *theRule     = theProblem->rule;
     femDiscrete    *theSpace    = theProblem->space;
     femGeo         *theGeometry = theProblem->geometry;
     femNodes       *theNodes    = theGeometry->theNodes;
     femMesh        *theMesh     = theGeometry->theElements;
     femMesh        *theEdges    = theGeometry->theEdges;
 
     double x[4], y[4], phi[4], dphidxsi[4] ,dphideta[4], dphidx[4], dphidy[4];
     int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];
     int nLocal = theMesh->nLocalNode;
     double a   = theProblem->A;
     double b   = theProblem->B;
     double c   = theProblem->C;      
     double rho = theProblem->rho;
     double g   = theProblem->g;
     double **A = theSystem->A;
     double *B  = theSystem->B;
     
     for (iElem = 0; iElem < theMesh->nElem; iElem++){
         for (j = 0; j < nLocal; j++){
             map[j]  = theMesh->elem[iElem * nLocal + j];
             mapX[j] = 2 * map[j];
             mapY[j] = 2 * map[j] + 1;
             x[j]    = theNodes->X[map[j]];
             y[j]    = theNodes->Y[map[j]];
         } 
         
         for (iInteg = 0; iInteg < theRule->n; iInteg++){   
             double xsi    = theRule->xsi[iInteg];
             double eta    = theRule->eta[iInteg];
             double weight = theRule->weight[iInteg];
 
             femDiscretePhi2(theSpace, xsi, eta, phi);
             femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
             
             double dxdxsi = 0.0; double dxdeta = 0.0;
             double dydxsi = 0.0; double dydeta = 0.0;
             for (i = 0; i < theSpace->n; i++){  
                 dxdxsi += x[i] * dphidxsi[i];       
                 dxdeta += x[i] * dphideta[i];   
                 dydxsi += y[i] * dphidxsi[i];   
                 dydeta += y[i] * dphideta[i];
             }
 
             double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
             
             for (i = 0; i < theSpace->n; i++){    
                 dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                 dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
             }
 
             double weightedJac = jac * weight;
 
             for (i = 0; i < theSpace->n; i++){ 
                 for (j = 0; j < theSpace->n; j++){
                     A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * weightedJac;
                     A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * weightedJac;
                     A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * weightedJac;
                     A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * weightedJac;
                 }
                 B[mapY[i]] -= phi[i] * g * rho * weightedJac;
             }
         }
     }
 }
 
void femElasticityAssembleNeumann(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->ruleEdge;
    femDiscrete    *theSpace = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[2],y[2],phi[2];
    int iBnd,iElem,iInteg,iEdge,i,j,d,map[2],mapU[2];
    int nLocal = 2;
    double *B  = theSystem->B;

    for(iBnd=0; iBnd < theProblem->nBoundaryConditions; iBnd++){
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        femDomain *theDomain = theCondition->domain;
        double value = theCondition->value;
        // Checking if the boundary condition is a Neumann condition
        if (type == NEUMANN_X || type == NEUMANN_Y){
            // Looping on the domain of the boundary condition
            for (iEdge = 0; iEdge < theDomain->nElem; iEdge++){
                // Getting the global element number
                iElem = theDomain->elem[iEdge];
                for (i = 0; i < nLocal; i++){
                    map[i] = theEdges->elem[iElem*nLocal + i];
                    mapU[i] = (type == NEUMANN_X) ? 2 * map[i] : 2 * map[i] + 1;
                    x[i] = theNodes->X[map[i]];
                    y[i] = theNodes->Y[map[i]];
                }
                // Computing the 1D jacobian
                double dx = x[1] - x[0];
                double dy = y[1] - y[0];
                double jacobian = sqrt(dx*dx + dy*dy)/2;
                // Integrating in 1D
                for (iInteg = 0; iInteg < theRule->n; iInteg++){
                    double xsi = theRule->xsi[iInteg];
                    double weight = theRule->weight[iInteg];
                    femDiscretePhi(theSpace, xsi, phi);
                    // Assembling Forces in B
                    for (i = 0; i < theSpace->n; i++){
                        B[mapU[i]] += phi[i] * value * jacobian * weight;
                    }
                }
            }
        }
    }
}

void femElasticityAssembleNeumannNormal(femProblem *theProblem){
    femFullSystem  *theSystem   = theProblem->system;
    femIntegration *theRule     = theProblem->ruleEdge;
    femDiscrete    *theSpace    = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes    = theGeometry->theNodes;
    femMesh        *theEdges    = theGeometry->theEdges;
    double x[2], y[2], phi[2];
    int iBnd, iElem, iInteg, iEdge, i, map[2];
    int nLocal = 2;
    double *B = theSystem->B;

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        femDomain *theDomain = theCondition->domain;

        if (type == NEUMANN_NORMAL) {
            for (iEdge = 0; iEdge < theDomain->nElem; iEdge++) {
                iElem = theDomain->elem[iEdge];

                for (i = 0; i < nLocal; i++) {
                    map[i] = theEdges->elem[iElem * nLocal + i];
                    x[i] = theNodes->X[map[i]];
                    y[i] = theNodes->Y[map[i]];
                }

                double dx = x[1] - x[0];
                double dy = y[1] - y[0];
                double norm = sqrt(dx*dx + dy*dy);
                double nx = -dy / norm;
                double ny =  dx / norm;
                double jacobian = norm / 2.0;

                for (iInteg = 0; iInteg < theRule->n; iInteg++) {
                    double xsi = theRule->xsi[iInteg];
                    double weight = theRule->weight[iInteg];
                    femDiscretePhi(theSpace, xsi, phi);

                    // Point d’intégration (x,y)
                    double xInteg = phi[0] * x[0] + phi[1] * x[1];
                    double yInteg = phi[0] * y[0] + phi[1] * y[1];

                    double value = (theCondition->profile != NULL)
                        ? theCondition->profile(xInteg, yInteg)
                        : theCondition->value;

                    double fx = value * nx;
                    double fy = value * ny;

                    for (i = 0; i < theSpace->n; i++) {
                        B[2*map[i]  ] += phi[i] * fx * jacobian * weight;
                        B[2*map[i]+1] += phi[i] * fy * jacobian * weight;
                    }
                }
            }
        }
    }
}
 
 
 double femMin(double *x, int n) 
 {
     double myMin = x[0];
     int i;
     for (i=1 ;i < n; i++) 
         myMin = fmin(myMin,x[i]);
     return myMin;
 }
 
 double femMax(double *x, int n) 
 {
     double myMax = x[0];
     int i;
     for (i=1 ;i < n; i++) 
         myMax = fmax(myMax,x[i]);
     return myMax;
 }
 
 double *femElasticitySolve(femProblem *theProblem, femSolverType FEM_GAUSS){
    femFullSystem *theSystem = theProblem->system;
    femFullSystemInit(theSystem);
    printf("Problem initiated\n");
    femElasticityAssembleElements(theProblem);
    printf("Elements assembled\n");
    femElasticityAssembleNeumann(theProblem);
    printf("Neumann assembled\n");
    femElasticityAssembleNeumannNormal(theProblem);
    printf("Neumann normal assembled\n");
    int size = theSystem->size;
    if (A_copy == NULL){
        A_copy = (double **) malloc(sizeof(double *) * size);
        for (int i = 0; i < size; i++) { A_copy[i] = (double *) malloc(sizeof(double) * size); }
    }
    if (B_copy == NULL) { B_copy = (double *) malloc(sizeof(double) * size); }

    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++) { A_copy[i][j] = theSystem->A[i][j]; }
        B_copy[i] = theSystem->B[i];
    }

    int *theConstrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < size; i++){
        if (theConstrainedNodes[i] != -1){
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem, i, value);
        }
    }

    printf("Solving the system\n");
    femFullSystemEliminate(theSystem, FEM_GAUSS);
    printf("System solved\n");
    // memcpy(theProblem->soluce, theSystem->B, theSystem->size * sizeof(double));
    // printf("Solution copied\n");
    // return theProblem->soluce;
    return theSystem->B;
 }
 
 double *femElasticityForces(femProblem *theProblem){
     double *residuals = theProblem->residuals;
     double *soluce    = theProblem->soluce;
     int size = theProblem->system->size;
 
     if (residuals == NULL) { residuals = (double *) malloc(sizeof(double) * size); }
 
     for (int i = 0; i < size; i++) { residuals[i] = 0.0; }
 
 
     for (int i = 0; i < size; i++){
         for (int j = 0; j < size; j++) { residuals[i] += A_copy[i][j] * soluce[j]; }
         residuals[i] -= B_copy[i];
     }
 
     for (int i = 0; i < size; i++) { free(A_copy[i]); A_copy[i] = NULL;}
     free(A_copy); free(B_copy);
     A_copy = NULL; B_copy = NULL;
 
     return residuals;
 }

 double foilProfile(double x, double y) {
    double amp = 5.0;
    double xShape = exp(-20.0 * (x - 0.1) * (x - 0.1));  // pic vers l'avant
    double yDir = (y < 0) ? -1.0 : 1.0;                  // vers le haut si dessous
    return amp * xShape * yDir;
}
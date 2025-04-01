#include "mesh.h"

femGeo theGeometry;

femGeo *geoGetGeometry(){ 
    return &theGeometry; 
}

double geoSizeDefault(double x, double y){ 
    return 1.0; 
}

double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data){
    return theGeometry.geoSize(x,y);    
}

double interpolationH(double x, double x0, double x1, double h0, double h) {
    double q0 = pow((x-x1)/(x0-x1),2);
    double q1 = pow((x-x0)/(x1-x0),2);
    double P0 = h0 + (x-x0)*(-2/(x0-x1))*h0;
    double P1 = h + (x-x1)*(-2/(x1-x0))*h;
    return q0*P0 + q1*P1;
}

double geoSize(double x, double y){
    femGeo *theGeometry = geoGetGeometry();
    double h = theGeometry->h;

    double x0 = theGeometry->xCircle1;
    double y0 = theGeometry->yCircle1;
    double r0 = theGeometry->rCircle1;
    double d0 = theGeometry->dCircle1;
    double h0 = theGeometry->hCircle1;

    double x1 = theGeometry->xCircle2;
    double y1 = theGeometry->yCircle2;
    double r1 = theGeometry->rCircle2;
    double d1 = theGeometry->dCircle2;
    double h1 = theGeometry->hCircle3;

    double x2 = theGeometry->xCircle3;
    double y2 = theGeometry->yCircle3;
    double r2 = theGeometry->rCircle3;
    double d2 = theGeometry->dCircle3;
    double h2 = theGeometry->hCircle3;

    double hfinal = h;
    double d = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) - r0;
    if (d < d0) {
        double a = (-2*h + 2*h0)/(d0*d0*d0);
        double b = (3*h  - 3*h0)/(d0*d0);
        double c = 0;
        hfinal = a*d*d*d + b*d*d + c*d + h0; }
        
    d = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1)) - r1;
    if (d < d1) {
        double a = (-2*h + 2*h1)/(d1*d1*d1);
        double b = (3*h  - 3*h1)/(d1*d1);
        double c = 0;
        hfinal = fmin(hfinal,a*d*d*d + b*d*d + c*d + h1); }
    
    d = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2)) - r2;
    if (d < d2) {
        double a = (-2*h + 2*h2)/(d2*d2*d2);
        double b = (3*h  - 3*h2)/(d2*d2);
        double c = 0;
        hfinal = fmin(hfinal,a*d*d*d + b*d*d + c*d + h2); }

    return hfinal;
}

int geoMeshGenerate() {
    int ierr;
    femGeo *theGeometry = geoGetGeometry();

    int N = theGeometry->N;
    double *joukowsky_x = theGeometry->joukowsky_x;
    double *joukowsky_y = theGeometry->joukowsky_y;

    double xCircle1 = theGeometry->xCircle1;
    double yCircle1 = theGeometry->yCircle1;
    double rCircle1 = theGeometry->rCircle1;

    double xCircle2 = theGeometry->xCircle2;
    double yCircle2 = theGeometry->yCircle2;
    double rCircle2 = theGeometry->rCircle2;

    double xCircle3 = theGeometry->xCircle3;
    double yCircle3 = theGeometry->yCircle3;
    double rCircle3 = theGeometry->rCircle3;

    // Wing
    int *points = malloc(N * sizeof(int));
    if (points == NULL) {printf("Error: Could not allocate memory for points.\n"); return 1;}
    for (int i = 0; i < N; i++) {
        points[i] = gmshModelOccAddPoint(joukowsky_x[i], joukowsky_y[i], 0, 0.1, -1, &ierr);
        if (ierr) {printf("Error: Could not add point %d.\n", i); return 1;}
    }
    points[N] = points[0];

    int splineID = gmshModelOccAddSpline(points, N+1, -1, NULL, 0, &ierr);
    if (ierr) {printf("Error: Could not add spline.\n"); return 1;}
    int loopID = gmshModelOccAddCurveLoop(&splineID, 1, -1, &ierr);
    if (ierr) {printf("Error: Could not add curve loop.\n"); return 1;}
    int surfaceID = gmshModelOccAddPlaneSurface(&loopID, 1, -1, &ierr);
    if (ierr) {printf("Error: Could not add plane surface.\n"); return 1;}

    int wing[] = {2, surfaceID};

    // Circles
    int circle1ID = gmshModelOccAddDisk(xCircle1,yCircle1,0.0,rCircle1,rCircle1,-1,NULL,0,NULL,0,&ierr); 
    if (ierr) {printf("Error: Could not add circle 1.\n"); return 1;}
    int circle2ID = gmshModelOccAddDisk(xCircle2,yCircle2,0.0,rCircle2,rCircle2,-1,NULL,0,NULL,0,&ierr);
    if (ierr) {printf("Error: Could not add circle 2.\n"); return 1;}
    int circle3ID = gmshModelOccAddDisk(xCircle3,yCircle3,0.0,rCircle3,rCircle3,-1,NULL,0,NULL,0,&ierr);
    if (ierr) {printf("Error: Could not add circle 3.\n"); return 1;}

    int circle1[] = {2, circle1ID};
    int circle2[] = {2, circle2ID};
    int circle3[] = {2, circle3ID};
    
    // Removing the circles from the wing
    gmshModelOccCut(wing,2,circle1,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    if (ierr) {printf("Error: Could not cut circle 1 from wing.\n"); return 1;}
    gmshModelOccCut(wing,2,circle2,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    if (ierr) {printf("Error: Could not cut circle 2 from wing.\n"); return 1;}
    gmshModelOccCut(wing,2,circle3,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    if (ierr) {printf("Error: Could not cut circle 3 from wing.\n"); return 1;}
    

    geoSetSizeCallback(geoSize);
    gmshModelOccSynchronize(&ierr);
    
    
    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",11,&ierr);  
        gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  
    }
    
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);
    }
    
    return 0;
}

void geoInitialize() {
    int ierr;
    theGeometry.geoSize = geoSizeDefault;
    gmshInitialize(0,NULL,1,0,&ierr);                         ErrorGmsh(ierr);
    gmshModelAdd("MyGeometry",&ierr);                         ErrorGmsh(ierr);
    gmshModelMeshSetSizeCallback(geoGmshSize,NULL,&ierr);     ErrorGmsh(ierr);
    theGeometry.theNodes = NULL;
    theGeometry.theElements = NULL;
    theGeometry.theEdges = NULL;
    theGeometry.nDomains = 0;
    theGeometry.theDomains = NULL;

    theGeometry.joukowsky_x = NULL;
    theGeometry.joukowsky_y = NULL;
}

void geoFinalize() {
    int ierr;
    
    if (theGeometry.theNodes) {
        free(theGeometry.theNodes->X);
        free(theGeometry.theNodes->Y);
        free(theGeometry.theNodes); }
    if (theGeometry.theElements) {
        free(theGeometry.theElements->elem);
        free(theGeometry.theElements); }
    if (theGeometry.theEdges) {
        free(theGeometry.theEdges->elem);
        free(theGeometry.theEdges); }
    for (int i=0; i < theGeometry.nDomains; i++) {
        free(theGeometry.theDomains[i]->elem);
        free(theGeometry.theDomains[i]);  }
    free(theGeometry.theDomains);

    if (theGeometry.joukowsky_x) free(theGeometry.joukowsky_x);
    if (theGeometry.joukowsky_y) free(theGeometry.joukowsky_y);

    gmshFinalize(&ierr); ErrorGmsh(ierr);
}

void geoSetSizeCallback(double (*geoSize)(double x, double y)) {
    theGeometry.geoSize = geoSize; 
}

void geoSetDomainName(int iDomain, char *name) {
    if (iDomain >= theGeometry.nDomains)  Error("Illegal domain number");
    if (geoGetDomain(name) != -1)         Error("Cannot use the same name for two domains");
    sprintf(theGeometry.theDomains[iDomain]->name,"%s",name);
} 

int geoGetDomain(char *name){
    int theIndex = -1;
    int nDomains = theGeometry.nDomains;
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        if (strncasecmp(name,theDomain->name,MAXNAME) == 0)
            theIndex = iDomain;  }
    return theIndex;        
}

void geoMeshImport() {
    int ierr;
    
    /* Importing nodes */
    
    size_t nNode,n,m,*node;
    double *xyz,*trash;
  //  gmshModelMeshRenumberNodes(&ierr);                        ErrorGmsh(ierr);
    gmshModelMeshGetNodes(&node,&nNode,&xyz,&n,
                         &trash,&m,-1,-1,0,0,&ierr);          ErrorGmsh(ierr);                         
    femNodes *theNodes = malloc(sizeof(femNodes));
    theNodes->nNodes = nNode;
    theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
    theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
    for (int i = 0; i < theNodes->nNodes; i++){
        theNodes->X[i] = xyz[3*node[i]-3];
        theNodes->Y[i] = xyz[3*node[i]-2]; }
    theGeometry.theNodes = theNodes;

    theNodes->number = malloc(sizeof(int)*theNodes->nNodes);
    for (int i = 0; i < theNodes->nNodes; i++) theNodes->number[i] = i;
    
    gmshFree(node);
    gmshFree(xyz);
    gmshFree(trash);
    printf("Geo     : Importing %d nodes \n",theGeometry.theNodes->nNodes);
       
    /* Importing elements */
    /* Pas super joli : a ameliorer pour eviter la triple copie */
        
    size_t nElem, *elem;
    gmshModelMeshGetElementsByType(1,&elem,&nElem,
                               &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
    femMesh *theEdges = malloc(sizeof(femMesh));
    theEdges->nLocalNode = 2;
    theEdges->nodes = theNodes;
    theEdges->nElem = nElem;  
    theEdges->elem = malloc(sizeof(int)*2*theEdges->nElem);
    for (int i = 0; i < theEdges->nElem; i++)
        for (int j = 0; j < theEdges->nLocalNode; j++)
            theEdges->elem[2*i+j] = node[2*i+j]-1;  
    theGeometry.theEdges = theEdges;
    int shiftEdges = elem[0];
    gmshFree(node);
    gmshFree(elem);
    printf("Geo     : Importing %d edges \n",theEdges->nElem);
  
    gmshModelMeshGetElementsByType(2,&elem,&nElem,
                               &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
    if (nElem != 0) {
      femMesh *theElements = malloc(sizeof(femMesh));
      theElements->nLocalNode = 3;
      theElements->nodes = theNodes;
      theElements->nElem = nElem;  
      theElements->elem = malloc(sizeof(int)*3*theElements->nElem);
      for (int i = 0; i < theElements->nElem; i++)
          for (int j = 0; j < theElements->nLocalNode; j++)
              theElements->elem[3*i+j] = node[3*i+j]-1;  
      theGeometry.theElements = theElements;
      gmshFree(node);
      gmshFree(elem);
      printf("Geo     : Importing %d triangles \n",theElements->nElem); }
    
    int nElemTriangles = nElem;
    gmshModelMeshGetElementsByType(3,&elem,&nElem,
                               &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
    if (nElem != 0 && nElemTriangles != 0)  
      Error("Cannot consider hybrid geometry with triangles and quads :-(");                       
                               
    if (nElem != 0) {
      femMesh *theElements = malloc(sizeof(femMesh));
      theElements->nLocalNode = 4;
      theElements->nodes = theNodes;
      theElements->nElem = nElem;  
      theElements->elem = malloc(sizeof(int)*4*theElements->nElem);
      for (int i = 0; i < theElements->nElem; i++)
          for (int j = 0; j < theElements->nLocalNode; j++)
              theElements->elem[4*i+j] = node[4*i+j]-1;  
      theGeometry.theElements = theElements;
      gmshFree(node);
      gmshFree(elem);
      printf("Geo     : Importing %d quads \n",theElements->nElem); }

    
    /* Importing 1D entities */
  
    int *dimTags;
    gmshModelGetEntities(&dimTags,&n,1,&ierr);        ErrorGmsh(ierr);
    theGeometry.nDomains = n/2;
    theGeometry.theDomains = malloc(sizeof(femDomain*)*n/2);
    printf("Geo     : Importing %d entities \n",theGeometry.nDomains);

    for (int i=0; i < n/2; i++) {
        int dim = dimTags[2*i+0];
        int tag = dimTags[2*i+1];
        femDomain *theDomain = malloc(sizeof(femDomain)); 
        theGeometry.theDomains[i] = theDomain;
        theDomain->mesh = theEdges;
        sprintf(theDomain->name, "Entity %d ",tag-1);
         
        int *elementType;
        size_t nElementType, **elementTags, *nElementTags, nnElementTags, **nodesTags, *nNodesTags, nnNodesTags; 
        gmshModelMeshGetElements(&elementType, &nElementType, &elementTags, &nElementTags, &nnElementTags, &nodesTags, &nNodesTags, &nnNodesTags, dim, tag, &ierr);
        theDomain->nElem = nElementTags[0];
        theDomain->elem = malloc(sizeof(int)*2*theDomain->nElem); 
        for (int j = 0; j < theDomain->nElem; j++) {
            theDomain->elem[j] = elementTags[0][j] - shiftEdges; }
        printf("Geo     : Entity %d : %d elements \n",i,theDomain->nElem);
        gmshFree(nElementTags);
        gmshFree(nNodesTags);
        gmshFree(elementTags);
        gmshFree(nodesTags);
        gmshFree(elementType); }
    gmshFree(dimTags);
 
    return;

}

void geoMeshPrint() {
    femNodes *theNodes = theGeometry.theNodes;
    if (theNodes != NULL) {
       printf("Number of nodes %d \n", theNodes->nNodes);
       for (int i = 0; i < theNodes->nNodes; i++) {
         printf("%6d : %6d : %14.7e %14.7e \n",i,theNodes->number[i],theNodes->X[i],theNodes->Y[i]); }}
    femMesh *theEdges = theGeometry.theEdges;
    if (theEdges != NULL) {
      printf("Number of edges %d \n", theEdges->nElem);
      int *elem = theEdges->elem;
      for (int i = 0; i < theEdges->nElem; i++) {
         printf("%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }}
    femMesh *theElements = theGeometry.theElements;
    if (theElements != NULL) {
      if (theElements->nLocalNode == 3) {
         printf("Number of triangles %d \n", theElements->nElem);
         int *elem = theElements->elem;
         for (int i = 0; i < theElements->nElem; i++) {
             printf("%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }}
      if (theElements->nLocalNode == 4) {
         printf("Number of quads %d \n", theElements->nElem);
         int *elem = theElements->elem;
         for (int i = 0; i < theElements->nElem; i++) {
             printf("%6d : %6d %6d %6d %6d\n",i,elem[4*i],elem[4*i+1],elem[4*i+2],elem[4*i+3]); }}}
    int nDomains = theGeometry.nDomains;
    printf("Number of domains %d\n", nDomains);
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
       femDomain *theDomain = theGeometry.theDomains[iDomain];
       printf("  Domain : %6d \n", iDomain);
       printf("  Name : %s\n", theDomain->name);
       printf("  Number of elements : %6d\n", theDomain->nElem);
       for (int i=0; i < theDomain->nElem; i++){
  //         if (i != theDomain->nElem  && (i % 10) != 0)  printf(" - ");
           printf("%6d",theDomain->elem[i]);
           if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) printf("\n"); }
       printf("\n"); }
   
   
 }
  
 void geoMeshWrite(const char *filename) {
    FILE* file = fopen(filename,"w");
  
    femNodes *theNodes = theGeometry.theNodes;
    fprintf(file, "Number of nodes %d \n", theNodes->nNodes);
    for (int i = 0; i < theNodes->nNodes; i++) {
       fprintf(file,"%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }
       
    femMesh *theEdges = theGeometry.theEdges;
    fprintf(file,"Number of edges %d \n", theEdges->nElem);
    int *elem = theEdges->elem;
    for (int i = 0; i < theEdges->nElem; i++) {
       fprintf(file,"%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }
       
    femMesh *theElements = theGeometry.theElements;
    if (theElements->nLocalNode == 3) {
       fprintf(file,"Number of triangles %d \n", theElements->nElem);
       elem = theElements->elem;
       for (int i = 0; i < theElements->nElem; i++) {
           fprintf(file,"%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }}
    if (theElements->nLocalNode == 4) {
       fprintf(file,"Number of quads %d \n", theElements->nElem);
       elem = theElements->elem;
       for (int i = 0; i < theElements->nElem; i++) {
           fprintf(file,"%6d : %6d %6d %6d %6d\n",i,elem[4*i],elem[4*i+1],elem[4*i+2],elem[4*i+3]); }}
      
    int nDomains = theGeometry.nDomains;
    fprintf(file,"Number of domains %d\n", nDomains);
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
       femDomain *theDomain = theGeometry.theDomains[iDomain];
       fprintf(file,"  Domain : %6d \n", iDomain);
       fprintf(file,"  Name : %s\n", theDomain->name);
       fprintf(file,"  Number of elements : %6d\n", theDomain->nElem);
       for (int i=0; i < theDomain->nElem; i++){
           fprintf(file,"%6d",theDomain->elem[i]);
           if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) fprintf(file,"\n"); }
       fprintf(file,"\n"); }
     
    fclose(file);
    geoMeshUnfuck();
 }
 
 void geoMeshRead(const char *filename) 
 {
    FILE* file = fopen(filename,"r");
    
    int trash, *elem;
    
    femNodes *theNodes = malloc(sizeof(femNodes));
    theGeometry.theNodes = theNodes;
    ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
    theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
    theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
    for (int i = 0; i < theNodes->nNodes; i++) {
        ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theNodes->X[i],&theNodes->Y[i]));} 
     
     theNodes->number = malloc(sizeof(int)*theNodes->nNodes);
     for (int i = 0; i < theNodes->nNodes; i++) theNodes->number[i] = i;
 
    femMesh *theEdges = malloc(sizeof(femMesh));
    theGeometry.theEdges = theEdges;
    theEdges->nLocalNode = 2;
    theEdges->nodes = theNodes;
    ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
    theEdges->elem = malloc(sizeof(int)*theEdges->nLocalNode*theEdges->nElem);
    for(int i=0; i < theEdges->nElem; ++i) {
         elem = theEdges->elem;
         ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash,&elem[2*i],&elem[2*i+1])); }
   
    femMesh *theElements = malloc(sizeof(femMesh));
    theGeometry.theElements = theElements;
    theElements->nLocalNode = 0;
    theElements->nodes = theNodes;
    char elementType[MAXNAME];  
    ErrorScan(fscanf(file, "Number of %s %d \n",elementType,&theElements->nElem));  
    if (strncasecmp(elementType,"triangles",MAXNAME) == 0) {
       theElements->nLocalNode = 3;
       theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
       for(int i=0; i < theElements->nElem; ++i) {
           elem = theElements->elem;
           ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", 
                     &trash,&elem[3*i],&elem[3*i+1],&elem[3*i+2])); }}
    if (strncasecmp(elementType,"quads",MAXNAME) == 0) {
       theElements->nLocalNode = 4;
       theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
       for(int i=0; i < theElements->nElem; ++i) {
           elem = theElements->elem;
           ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", 
                     &trash,&elem[4*i],&elem[4*i+1],&elem[4*i+2],&elem[4*i+3])); }}
            
    ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
    int nDomains = theGeometry.nDomains;
    theGeometry.theDomains = malloc(sizeof(femDomain*)*nDomains);
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
       femDomain *theDomain = malloc(sizeof(femDomain)); 
       theGeometry.theDomains[iDomain] = theDomain;
       theDomain->mesh = theEdges; 
       ErrorScan(fscanf(file,"  Domain : %6d \n", &trash));
       ErrorScan(fscanf(file,"  Name : %[^\n]s \n", (char*)&theDomain->name));
       ErrorScan(fscanf(file,"  Number of elements : %6d\n", &theDomain->nElem));
       theDomain->elem = malloc(sizeof(int)*2*theDomain->nElem); 
       for (int i=0; i < theDomain->nElem; i++){
           ErrorScan(fscanf(file,"%6d",&theDomain->elem[i]));
           if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) ErrorScan(fscanf(file,"\n")); }}
     
    fclose(file);
 }

 void geoMeshUnfuck() {
    int bettercallsaul = system("fixmesh.py");
    if (bettercallsaul != 0) {
        printf("Error: fixmesh.py failed with error code %d\n", bettercallsaul);
        return;
    }
    if (bettercallsaul == 0) {
        printf("fixmesh.py executed successfully.\n");
    }
 }
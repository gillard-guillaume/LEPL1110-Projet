# include "../fem/fem.h"

double *toCompare;
int compareNodePosition(const void *a, const void *b){
    int i = *(int*)a;
    int j = *(int*)b;
    if (toCompare[i] < toCompare[j]) return 1;
    if (toCompare[i] > toCompare[j]) return -1;
    return 0;
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType){
    int i;
    int nNodes = theMesh->nodes->nNodes;
    int *temp = (int *) malloc(nNodes*sizeof(int));
    for (i = 0; i < nNodes; i++) {temp[i] = i;}
    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[i] = i;
            break;

        case FEM_XNUM :
            toCompare = theMesh->nodes->X;
            qsort(temp, nNodes, sizeof(int), compareNodePosition);
            break;

        case FEM_YNUM : 
            toCompare = theMesh->nodes->Y;
            qsort(temp, nNodes, sizeof(int), compareNodePosition);
            break;   

        default : Error("Unexpected renumbering option"); 
    }

    for (i = 0; i < nNodes; i++) 
        theMesh->nodes->number[temp[i]] = i;
    free(temp);
}

int femMeshComputeBand(femMesh *theMesh){
    int myBand = 0;
    int nElem = theMesh->nElem;
    int nLocalNode = theMesh->nLocalNode;
    int min, max;
    int map[nLocalNode];

    for (int iElem=0; iElem < nElem; iElem++){
        for (int j=0; j < nLocalNode; j++){map[j] = theMesh->elem[iElem*nLocalNode+j];}
        min = map[0];
        max = map[0];
        for (int j=1; j < nLocalNode; j++){
            min = MIN(min,map[j]);
            max = MAX(max,map[j]);
        }
        myBand = MAX(myBand,max-min+1);
    }
    return(myBand);
}
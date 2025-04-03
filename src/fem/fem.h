
/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <time.h>
#include "gmshc.h"


#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define ErrorGmsh(a)   femErrorGmsh(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1
#define MAXNAME 256
#define EPS 1e-15
#define MAX_ITER 1000

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef MAX
    #define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
    #define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;
typedef enum {FEM_TRIANGLE,FEM_QUAD,FEM_EDGE} femElementType;
typedef enum {DIRICHLET_X,DIRICHLET_Y,NEUMANN_X,NEUMANN_Y, NEUMANN_NORMAL} femBoundaryType;
typedef enum {PLANAR_STRESS,PLANAR_STRAIN,AXISYM} femElasticCase;
typedef enum {FEM_CG, FEM_CHOV, FEM_GAUSS} femSolverType;


typedef struct {
    int nNodes;
    double *X;
    double *Y;
    int *number;
} femNodes;

typedef struct {
    int nLocalNode;
    int nElem;
    int *elem;
    femNodes *nodes;
} femMesh;

typedef struct {
    femMesh *mesh;
    int nElem;
    int *elem;
    char name[MAXNAME];
} femDomain;

typedef struct {
    double R, muX, muY, N;
    double *joukowsky_x;
    double *joukowsky_y;
    double xCircle1, yCircle1, rCircle1, hCircle1, dCircle1;
    double xCircle2, yCircle2, rCircle2, hCircle2, dCircle2;
    double xCircle3, yCircle3, rCircle3, hCircle3, dCircle3;
    double h;
    femElementType elementType;
    double (*geoSize)(double x, double y);
    femNodes *theNodes;
    femMesh  *theElements;
    femMesh  *theEdges;
    int nDomains;
    femDomain **theDomains;
} femGeo;

typedef struct {
    int n;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
    void (*x)(double *xsi);
    void (*phi)(double xsi, double *phi);
    void (*dphidx)(double xsi, double *dphidxsi);
} femDiscrete;
    
typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;


typedef struct {
    femDomain* domain;
    femBoundaryType type; 
    double value;
} femBoundaryCondition;


typedef struct {
    double E,nu,rho,g;
    double A,B,C;
    int planarStrainStress;
    int nBoundaryConditions;
    femBoundaryCondition **conditions;  
    int *constrainedNodes; 
    double *soluce;
    double *residuals;
    femGeo *geometry;
    femDiscrete *space;
    femIntegration *rule;
    femDiscrete *spaceEdge;
    femIntegration *ruleEdge;
    femFullSystem *system;
} femProblem;

int geoGetDomain(char *name);
void geoSetDomainName(int iDomain, char *name);

femProblem*         femElasticityCreate(femGeo* theGeometry, 
                                      double E, double nu, double rho, double g, 
                                      femElasticCase iCase, femRenumType renumType);
void                femElasticityFree(femProblem *theProblem);
void                femElasticityPrint(femProblem *theProblem);
void                femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value);
void                femElasticityAssembleElements(femProblem *theProblem);
void                femElasticityAssembleNeumann(femProblem *theProblem);
double*             femElasticitySolve(femProblem *theProblem, femSolverType solverType);
double*             femElasticityForces(femProblem *theProblem);
double              femElasticityIntegrate(femProblem *theProblem, double (*f));

femIntegration*     femIntegrationCreate(int n, femElementType type);
void                femIntegrationFree(femIntegration *theRule);

femDiscrete*        femDiscreteCreate(int n, femElementType type);
void                femDiscreteFree(femDiscrete* mySpace);
void                femDiscretePrint(femDiscrete* mySpace);
void                femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);
void                femDiscreteXsi(femDiscrete* mySpace, double *xsi);
void                femDiscretePhi(femDiscrete* mySpace, double xsi, double *phi);
void                femDiscreteDphi(femDiscrete* mySpace, double xsi, double *dphidxsi);


femFullSystem*      femFullSystemCreate(int size);
void                femFullSystemFree(femFullSystem* mySystem);
void                femFullSystemPrint(femFullSystem* mySystem);
void                femFullSystemInit(femFullSystem* mySystem);
void                femFullSystemAlloc(femFullSystem* mySystem, int size);
double*             femFullSystemEliminate(femFullSystem* mySystem, femSolverType solverType);
void                femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);

double              femMin(double *x, int n);
double              femMax(double *x, int n);
void                femError(char *text, int line, char *file);
void                femErrorScan(int test, int line, char *file);
void                femErrorGmsh(int test, int line, char *file);
void                femWarning(char *text, int line, char *file);
void                joukovskyError(char *text, int line, char *file);

double complex      joukowskyTransform(double complex zeta);
double complex      inverseJoukowskyTransform(double complex z);
int                 joukowsky(femGeo *theGeometry);

int                 femMeshComputeBand(femMesh *theMesh);
void                femMeshRenumber(femMesh *theMesh, femRenumType renumType);

void                conjugateGradient(double **A, double *b, int n);
void                cholevsky(double **A, double *B, int size);
void                gauss(double **A, double *B, int size);



#endif
# include "fem.h"

double norm(double *v, int n){
    double sum = 0;
    for (int i=0; i<n; i++){
        sum += v[i]*v[i];
    }
    return sqrt(sum);
}

void conjugateGradient(double **A, double *b, int n){
    // Initialisation des variables
    double *r = (double *)malloc(n*sizeof(double));
    double *p = (double *)malloc(n*sizeof(double));

    double alpha;
    double beta;
    double rtr;
    double new_rtr;
    double *Ap = (double *)malloc(n*sizeof(double));
    double pAp;

    // On stock la solution dans b
    for (int i = 0; i<n; i++){
        r[i] = b[i];
        p[i] = b[i];
        b[i] = 0;
    }

    // Gradient conjugué
    rtr = 0;
    for (int i=0; i<n; i++) rtr += r[i]*r[i];

    for (int i=0; i<MAX_ITER; i++){
        // A p
        for (int j=0; j<n; j++){
            Ap[j] = 0;
            for (int k=0; k<n; k++){
                Ap[j] += A[j][k]*p[k];}}

        // p^T A p
        pAp = 0;
        for (int j=0; j<n; j++) pAp += p[j]*Ap[j];
        // Longueur du pas
        alpha = rtr/pAp;

        // Solution approximative
        for (int j=0; j<n; j++) b[j] += alpha*p[j];

        // Résidu
        for (int j=0; j<n; j++) r[j] -= alpha*Ap[j];
        
        // Si le résidu est suffisament petit, on s'arrête
        if (norm(r, n) < EPS) break;

        // Amélioration du pas
        new_rtr = 0;
        for (int j=0; j<n; j++) new_rtr += r[j]*r[j];
        beta = new_rtr/rtr;
        rtr = new_rtr;

        // Nouvelle direction
        for (int j=0; j<n; j++) p[j] = r[j] + beta*p[j];
        
        if (norm(r, n) < EPS){
            printf("Convergence atteinte en %d itérations\n", i);
            break;
        }
    }
    free(r);
    free(p);
    free(Ap);
    return;
}

void backSubstitution(double **A, double *B, int size){
    int factor;
    for (int i = size-1; i >= 0 ; i--) {
        factor = 0;
        for (int j = i+1 ; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
}

void cholevsky(double **A, double *B, int size){
    int factor;
    for (int k=0; k<size; k++){
        for (int i = k+1; i<size; i++){
            double factor = A[k][i]/A[k][k];
            for (int j=i; j<size; j++){
                A[i][j] -= factor*A[k][j];
            }
            B[i] -= factor*B[k];
        }
    }
    backSubstitution(A, B, size);
}

void gauss(double **A, double *B, int size){
    int factor;
    for (int k=0; k < size; k++) {
        for (int i = k+1 ; i <  size; i++) {
            factor = A[i][k] / A[k][k];
            for (int j = k+1 ; j < size; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; 
        }
    }

    backSubstitution(A, B, size);
}
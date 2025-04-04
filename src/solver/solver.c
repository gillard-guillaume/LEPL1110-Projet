# include "../fem/fem.h"

double norm(double *v, int n){
    double sum = 0;
    for (int i=0; i<n; i++){
        sum += v[i]*v[i];
    }
    return sqrt(sum);
}


void conjugateGradient(double **A, double *b, double *x, int n) {
    double *r = malloc(n * sizeof(double));
    double *d = malloc(n * sizeof(double));
    double *Ad = malloc(n * sizeof(double));
    
    int i, j, k;
    double alpha, beta, rtr, new_rtr, dad;

    // Initialisation : r = b - Ax, d = r, x initialisé à 0
    for (i = 0; i < n; i++) {
        r[i] = b[i];
        for (j = 0; j < n; j++) {
            r[i] -= A[i][j] * x[j]; // Ax initialisé à 0 donc r = b
        }
        d[i] = r[i]; // d = r
    }

    rtr = 0;
    for (i = 0; i < n; i++) rtr += r[i] * r[i];

    for (k = 0; k < MAX_ITER; k++) {
        // Calcul de A * d
        for (i = 0; i < n; i++) {
            Ad[i] = 0;
            for (j = 0; j < n; j++) {
                Ad[i] += A[i][j] * d[j];
            }
        }

        // Calcul de d^T * Ad
        dad = 0;
        for (i = 0; i < n; i++) dad += d[i] * Ad[i];

        // Vérification de la division par zéro
        if (fabs(dad) < 1e-10) {
            printf("Erreur : d^T * A * d est trop proche de zéro. Arrêt.\n");
            break;
        }

        // Calcul du pas alpha
        alpha = rtr / dad;

        // Mise à jour de x : x = x + alpha * d
        for (i = 0; i < n; i++) x[i] += alpha * d[i];

        // Mise à jour de r : r = r - alpha * Ad
        for (i = 0; i < n; i++) r[i] -= alpha * Ad[i];

        // Vérification de la convergence
        if (norm(r, n) < EPS) {
            printf("Convergence atteinte en %d itérations\n", k + 1);
            break;
        }

        // Mise à jour du coefficient beta
        new_rtr = 0;
        for (i = 0; i < n; i++) new_rtr += r[i] * r[i];
        beta = new_rtr / rtr;
        rtr = new_rtr;

        // Mise à jour de la direction de recherche : d = r + beta * d
        for (i = 0; i < n; i++) d[i] = r[i] + beta * d[i];
    }

    // Libération de la mémoire
    free(r);
    free(d);
    free(Ad);
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
    double factor;
    int i, j, k;
    for (k=0; k < size; k++) {
        for (i = k+1 ; i <  size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; 
        }
    }
    backSubstitution(A, B, size);
}
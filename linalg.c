//
// Created by Pallas Dimitrios on 1/4/2022.
//

#include <stdlib.h>
#include "math.h"
#include <stdio.h>

void print_vector(int N , double X[N]){
    for(int i = 0 ; i < N ; i++){
        printf("%lf\n",X[i]);
    }
}

void print_matrix(int N , double X[N][N]){
    for(int i = 0 ; i < N ; i++) {
        for (int j = 0; j < N; j++) {
            printf(" %lf", X[i][j]);
        }
        printf("\n");
    }
}

void scalar_times_vector(int N , double scalar , double X[N]){
    for(int i = 0 ; i < N ; i++){
        X[i] = X[i] * scalar;
    }
}
void copy_matrix(int rows , int cols , double original[][cols] , double copy[][cols]){
    register int i ; register int j;
    for( i = 0 ; i < rows ; i++){
        for( j = 0 ; j < cols ; j++){
            copy[i][j] = original[i][j];
        }
    }
}

void copy_vector(int n , double original[n] , double copy[n]){
    register int i;
    for( i = 0 ; i < n ; i++){
        copy[i] = original[i];
    }
}
//Gauss elimination and back-subsitution
void solve_system(int n , double A_matrix[][n], double B_vector[n]  , double BM[n]){
    double AM[n][n];
    copy_matrix(n,n,A_matrix,AM);
    copy_vector(n,B_vector,BM);
    double scaler ;
    register int fd;
    register int i ;
    register int k ;
    register int j ;
    register double current_scaler;
    for( fd = 0 ; fd < n ; fd++) {
        scaler = 1 / (AM[fd][fd]);
        for ( j = 0; j < n; j++) {
            AM[fd][j] *= scaler;
        }
        BM[fd] *= scaler;
        for ( i = 0; i < fd; i++) {
            current_scaler = AM[i][fd];
            for ( k = 0; k < n; k++) {
                AM[i][k] = AM[i][k] - current_scaler * AM[fd][k];
            }
            BM[i] = BM[i] - current_scaler * BM[fd];
        }
        for ( i = fd + 1; i < n; i++) {
            current_scaler = AM[i][fd];
            for ( k = 0; k < n; k++) {
                AM[i][k] -= current_scaler * AM[fd][k];
            }
            BM[i] -= current_scaler * BM[fd];
        }
    } // Now the matrix AM is morphed into a nxn identity matrix
    // while BM contains a vector form solution of the system
}
double L2(int N , double X[N]){
    double l2norm = 0 ;
    for(int i = 0 ; i < N ; i++){
        l2norm += pow((X[i]),2);
    }
    return sqrt(l2norm);
}
double L2DIFF(int N , double X1[N] , double X2[N]){
    double l2norm = 0 ;
    for(int i = 0 ; i < N ; i++){
        l2norm += pow((X1[i]-X2[i]),2);
    }
    return sqrt(l2norm);
}


//Only use for systems where the matrix is strictly diagonally dominant or symmetric and postive definite
void gauss_seidel(int N , double A[N][N] , double B[N] , double X[N]){

    double X_new[N];
    double error[N];
    double L2NORM;
    for(int iter = 0; iter < 10000 ; iter++){
    for(int i = 0; i < N ; i++){
        X_new[i] = B[i] ;
        for(int j = 0; j < N ; j++) {
            if (i != j) {
                X_new[i] -= X[i] / A[i][j];
            }
        }
        X_new[i] = 0.5*X_new[i]/A[i][i] - 0.5*X[i];
        error[i] = X_new[i] - X[i];
    }L2NORM = L2(N, error);
    if(L2NORM<1e-6){
        break;
    }
        copy_vector(N,X_new,X);
    }
printf("no convergence");
}





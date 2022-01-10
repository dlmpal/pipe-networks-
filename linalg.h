//
// Created by Pallas Dimitrios on 1/4/2022.
//

#ifndef PIPENETWORK_LINALG_H
#define PIPENETWORK_LINALG_H
void print_matrix(int N , double X[N][N]);
void print_vector(int N , double X[N]);
void scalar_times_vector(int N , double scalar , double X[N]);
void solve_system(int n , double A_matrix[][n], double B_vector[n]  , double BM[n]);
void gauss_seidel(int N , double A[N][N] , double B[N] , double X[N]);
void copy_vector(int n , double original[n] , double copy[n]);
double L2DIFF(int N , double X1[N] , double X2[N]);
double L2(int N , double X[N]);
#endif //PIPENETWORK_LINALG_H

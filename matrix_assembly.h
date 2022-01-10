//
// Created by Pallas Dimitrios on 1/4/2022.
//

#ifndef PIPENETWORK_MATRIX_ASSEMBLY_H
#define PIPENETWORK_MATRIX_ASSEMBLY_H
void create_rhs(int n_NODES,int n_ELE,int ELE[n_ELE][2] ,int BC ,
                double H[n_NODES] , double K[n_ELE] , double F[n_NODES],double rhs[n_NODES]);

void create_Jacobian(int n_NODES , int n_ELE , int ELE[n_ELE][2] , int BC ,
                     double H[n_NODES] ,double K[n_ELE] , double J[n_NODES][n_NODES]);

void newton_solve(int n_NODES , int n_ELE , int ELE[n_ELE][2] , int BC, double F[n_NODES],
                  double L[n_ELE] , double D[n_ELE] , double Z[n_ELE],
                  double K[n_ELE] , double U[n_ELE] , double Q[n_ELE],
                  double H[n_NODES] , int max_iters , double tol  );

void create_K(int init , int n_ELE , double U[n_ELE] , double D[n_ELE] , double L[n_ELE] , double Z[n_ELE] , double K[n_ELE]);
void update_diameters(int n_ELE , int ELE[][2] , double U[n_ELE] , double Q[n_ELE] , double D[n_ELE]);
void calc_static(int n_NODES , int n_ELE , int ELE[][2] , double H[n_NODES] ,double h[n_NODES],double U[n_ELE] , double H_static[n_NODES]);
double cost_function(int n_ELE , double D[n_ELE] , double U[n_ELE]);

#endif //PIPENETWORK_MATRIX_ASSEMBLY_H

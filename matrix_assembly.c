//
// Created by Pallas Dimitrios on 1/4/2022.
//
#include <stdlib.h>
#include "math.h"
#include <stdio.h>
#include <string.h>
#include "linalg.h"

// defines the flow direction
double sign(double P1,double P2){
    if(P1>P2){
        return -1;
    }
    else{
        return 1;
    }
}
// defines the signs of the elements of the rhs vector
double sign2(int i , int N1 , int N2 , double *H , int BC){
    if(i == BC){
        return 0;
    }
    if(i != N1 && i != N2){
        return 0;
    }
    if(i == N1){
        return sign(H[N1],H[N2]);
    }
    if(i == N2){
        return sign(H[N2],H[N1]);
    }
    else{
        printf("Error in Connectivity at node %d",i);
        exit(-1);
    }
}
// defines the signs of the elements of the Jacobian
double sign3(int i, int j, int N1, int N2) {
    if (i != N1 && i != N2) {
        return 0;
    } else {
        if (i == j) {
            return -1;
        } else {
            if (i == N1 && j == N2) {
                return 1;
            }

            if (i == N2 && j == N1) {
                return 1;
            } else {
                return 0;
            }
        }
    }
}
// creates the right hand side vector (J@X = F, F = rhs)
void create_rhs(int n_NODES,int n_ELE,int ELE[n_ELE][2] ,int BC ,
                double H[n_NODES] , double K[n_ELE] , double F[n_NODES],double rhs[n_NODES]){
    //double *rhs;
    //rhs = malloc(sizeof(double)*n_NODES);
    for(int i = 0 ; i < n_NODES ;i++){
        rhs[i] = F[i];
        for(int ele = 0 ; ele < n_ELE ; ele++ ){
            rhs[i] += sqrt(fabs(H[ELE[ele][0]]-H[ELE[ele][1]])/K[ele]) * sign2(i , ELE[ele][0],ELE[ele][1],H,BC);
        }
    }
}
// creates the Jacobian matrix
void create_Jacobian(int n_NODES , int n_ELE , int ELE[n_ELE][2] , int BC , double H[n_NODES] ,double K[n_ELE] , double J[n_NODES][n_NODES]){
    for(int i = 0 ; i < n_NODES ; i++){
        for(int j = 0 ; j < n_NODES ; j++){
            if(i == BC){
                J[i][j] = (double) 1.0 * (i==j);
                continue;
            }
            for(int ele = 0 ; ele < n_ELE ; ele++){
                J[i][j] += 1/(2*sqrt(fabs(H[ELE[ele][0]]-H[ELE[ele][1]])*K[ele])) * sign3(i,j,ELE[ele][0],ELE[ele][1]);
            }
        }
    }
}
// calculates loss factor using Jain's formula
double jain(double u , double D){
    double epsilon = 5e-5;
    double nu = 13.57 * 1e-6;
    double Re = u*D/nu;
    return pow(1.14 - 2* log10(epsilon/D + 21.25/(pow(Re,0.9))),-2);
}
// the same only for Re>>1
double jain_0(double D){
    double epsilon = 5e-5;
    return pow(1.14-2*log10(epsilon/D),-2);
}

//creates the K matrix given the required args
void create_K(int init , int n_ELE , double U[n_ELE] , double D[n_ELE] , double L[n_ELE] , double Z[n_ELE] , double K[n_ELE]){
    for(int ele = 0 ; ele < n_ELE ; ele++){
        double zeta_lin;
        if(init){
            zeta_lin = jain_0(D[ele]);
        }
        else{
            zeta_lin = jain(U[ele],D[ele]);
        }
        K[ele] = 8/(9.81 * M_PI*M_PI*pow(D[ele],4))*(Z[ele] +zeta_lin*L[ele]/D[ele]);
    }
}
//updates all pipe velocities after a newton iteration
void update_U(int n_NODES , int n_ELE , int ELE[n_ELE][2] , double H[n_NODES],
              double D[n_ELE],double L[n_ELE] , double Z[n_ELE] , double K[n_ELE] , double Q[n_ELE] , double U[n_ELE] ){
    for(int ele = 0; ele < n_ELE ; ele++){
        Q[ele] = sqrt(fabs(H[ELE[ele][1]] - H[ELE[ele][0]]) / K[ele]);
        U[ele] = Q[ele] / (M_PI * pow(D[ele],2)/4);
    }
    create_K(0 , n_ELE , U , D , L , Z , K );
}


//defines and solves the non-linear system modelling the pipe network
void newton_solve(int n_NODES , int n_ELE , int ELE[n_ELE][2] , int BC, double F[n_NODES],
                 double L[n_ELE] , double D[n_ELE] , double Z[n_ELE],
                 double K[n_ELE] , double U[n_ELE] , double Q[n_ELE],
                 double H[n_NODES] , int max_iters , double tol  ){
    double corr_log[100];
    memset(corr_log,0,100*sizeof (double));
    create_K(1,n_ELE,U,D,L,Z,K);
    double dH[n_NODES];

    for(int iter = 0; iter <max_iters ; iter++){
        double L2NORM;
        double rhs[n_NODES];
        double jacobian[n_NODES][n_NODES];
        memset(jacobian,0,sizeof(double)*n_NODES*n_NODES);
        create_rhs(n_NODES,n_ELE,ELE,BC,H,K,F,rhs);
        create_Jacobian(n_NODES,n_ELE,ELE,BC,H,K,jacobian);
        solve_system(n_NODES,jacobian,rhs,dH);
        for(int i = 0; i < n_NODES ; i++){
            H[i] -= dH[i];
        }
        update_U(n_NODES,n_ELE,ELE,H,D,L,Z,K,Q,U);
        L2NORM = L2(n_NODES,dH);
        corr_log[iter] = L2NORM;
        if(L2NORM<tol){
            //printf("Newton Iterations : %d\n",iter+1);
            //print_vector(iter,corr_log);
            break;
        }
    }
}

// finds the max velocity of all elements conencted to a node
void max_neighbour_vel(int n_Nodes , int n_ELE , int ELE[][2] , double U[n_ELE] , double max_U[n_Nodes]){

    for(int i = 0 ; i < n_Nodes ; i++){
        double max_vel = 0;
        for(int ele = 0 ; ele < n_ELE ; ele++){
            if(ELE[ele][0]==i || ELE[ele][1]==i){
               if(max_vel < U[ele]) {
                   max_vel = U[ele];
               }
            }
        }
        max_U[i] = max_vel;
    }
}

//calculates static pressure
void calc_static(int n_NODES , int n_ELE , int ELE[][2] , double H[n_NODES] ,
                 double h[n_NODES],double U[n_ELE] , double H_static[n_NODES]){
    double rho_air = 1.225;
    double rho_gas = 0.79;
    H_static[0] = H[0];
    double max_U[n_NODES];
    max_neighbour_vel(n_NODES,n_ELE,ELE,U,max_U);
    for(int i = 1; i < n_NODES;i++){
        H_static[i] = H[i] - 0.5*rho_gas*max_U[i]*max_U[i]-(rho_air-rho_gas)*9.81*h[i];
    }
}

// finds the index of the closest standardized diameter
int match_diameter(double diameter){
    int N = 15;
    double std_diameters[] = {21.3-2.6,26.9-2.6,33.7-2.6,42.4-2.6,48.3-2.6,60.3-2.9,76.1-2.9,88.9-3.2,114.3-3.6,139.7-4,168.3-4.5,219.1-5.9
            ,273.0-6.3,323.9-7.1,355.6-7.1};
    double min_error , cur_error;
    min_error = 1e5;
    int count = 0;
    for(int i = 0 ; i < N ; i++){
        cur_error = fabs(diameter - std_diameters[i]*1e-3) / (std_diameters[i]*1e-3);
        if(cur_error < min_error){
            min_error = cur_error;
            count = i;
        }
    }
    return count;
}
//Used to upadte pipe diameters after a call to newton_solve
void update_diameters(int n_ELE , int ELE[][2] , double U[n_ELE] , double Q[n_ELE] , double D[n_ELE]){
    int N = 15;
    double std_diameters[] = {21.3-2.6,26.9-2.6,33.7-2.6,42.4-2.6,48.3-2.6,60.3-2.9,76.1-2.9,88.9-3.2,114.3-3.6,139.7-4,168.3-4.5,219.1-5.9
            ,273.0-6.3,323.9-7.1,355.6-7.1};
    for(int ele = 0 ; ele < n_ELE ; ele++){
        while(1){
            if(U[ele] > 3){
                int index = match_diameter(D[ele]);
                D[ele] = std_diameters[index+1]*1e-3;
                U[ele] = Q[ele]/(M_PI*(D[ele]*D[ele])/4);
                continue;
            }
            if(U[ele] < 1){
                int index = match_diameter(D[ele]);
                if(index==0){
                    break;
                }
                D[ele] = std_diameters[index-1]*1e-3;
                U[ele] = Q[ele]/(M_PI*(D[ele]*D[ele])/4);
                continue;
            }
            break;
        }
    }
}
//Not used-Can be used for optimization through adjoint/FD
double cost_function(int n_ELE , double D[n_ELE] , double U[n_ELE]){
    double cost = 0;
    for(int ele = 0 ; ele < n_ELE ; ele++){
    if(U[ele] > 6 || U[ele] < 0.5){
        cost += 1000;
    }
    cost += D[ele]*100;
    }
    return cost;
}


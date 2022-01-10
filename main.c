#include <stdio.h>
#include "matrix_assembly.h"
#include <string.h>
#include "linalg.h"
#include <time.h>

int main() {
    //calculate total program runtime
    clock_t time_start , timestop;
    time_start = clock();

    int I=2; // student First and Last name
    // some boilerplate-defining the network
    int n_NODES = 20;
    int n_ELE = 23;

    // element/pipe 1-2 is conncted to nodes 1 and 2
    int ELE[][2] = {
            {1,2},
            {2,3},
            {3,4},
            {4,5},
            {5,6},
            {6,7},
            {7,8},
            {8,9},
            {9,10},
            {3,10},
            {7,11},
            {3,11},
            {3,12},
            {12,13},
            {13,14},
            {14,15},
            {15,16},
            {16,17},
            {17,18},
            {18,19},
            {12,19},
            {12,20},
            {16,20},
    };
    //node 1 is node 0 for the computer!
    for(int i = 0 ; i < n_ELE ; i++){
        ELE[i][0] -=1;
        ELE[i][1] -=1;
    }
    //more boilerplate
    double D0 = (60-3.2)*1e-3;
    double Q_total = (2*(1.5+1.5+3+2)*0.234 + 2*(2+3)*0.325 + 2*15*0.8)/3600 * I;
    double F[] = {0,0,0,-1.5*0.234 , -1.5*0.234 , -2*0.325 , 0 , -3 * 0.325 , -3 * 0.234 , -2 * 0.234 , -15 * 0.8,
            0,-1.5*0.234 , -1.5*0.234 , -2*0.325 , 0 , -3 * 0.325 , -3 * 0.234 , -2 * 0.234 , -15 * 0.8};
    scalar_times_vector(n_NODES,I/3600.0,F);
    //initial choice of pressures
    double H[] = {2500,
            2490.0,
            2485.0,
            2480.0,
            2475.0,
            2470.0,
            2465.0,
            2460.0,
            2455.0,
            2450.0,
            2445.0,
            2440.0,
            2435.0,
            2430.0,
            2425.0,
            2420.0,
            2415.0,
            2410.0,
            2405.0,
            2400.0};
    //the altitude of each node
    double h[] = {0 , 0 , 3.5 , 3.5 , 3.5 , 3.5 , 3.5 , 3.5 , 3.5 , 3.5 ,
                  3.5 , 6.5 , 6.5 , 6.5 , 6.5, 6.5,6.5,6.5,6.5,6.5};
    scalar_times_vector(n_NODES,1/(9.81*0.79),H);
    //point losses
    double Z[] =
            {2.5, 3, 4.5, 2.8, 2, 1.3,
             1.3, 2, 2.8, 4.5, 1.5, 4, 3.8,
             4.5, 2.8, 2, 1.3,
             1.3, 2, 2.8, 4.5, 4, 1.5};
    //length of each pipe
    double L[] = {5, 3.5, 16, 5, 5, 16, 16,
                   5, 5, 16, 5, 5, 3, 16, 5,
                   5, 16, 16, 5, 5, 16, 5, 5};
    double D[n_ELE];
    for(int i = 0 ; i < n_ELE;i++){
        D[i] = D0;
    }
    double U[n_ELE],Q[n_ELE],K[n_ELE];
    memset(U,0,sizeof(double)*n_ELE);
    memset(Q,0,sizeof(double)*n_ELE);

    int system_iterations = 5; // how many times to call the newton_solve function
    for(int iter = 0 ; iter < system_iterations ; iter++) {
        newton_solve(n_NODES,n_ELE,ELE,0,F,L,D,Z,K,U,Q,H,1000,1e-6);
        if(iter==system_iterations-1){
            scalar_times_vector(n_NODES,0.79*9.81,H);
           // print_vector(n_NODES,H);
            scalar_times_vector(n_NODES,1/0.79/9.81,H);

        }
        if(iter == system_iterations -1){
            break;
        }
        update_diameters(n_ELE,ELE,U,Q,D);
    }
    printf("Final Diameters:\n");
    print_vector(n_ELE,D);
    printf("Final Velocities:\n");
    print_vector(n_ELE,U);

    //Network Test at fluctuations of daily demand
    char fluctuations_filepath[] = "C:/Users/dlmpa/Numerical_data/5th_semester/AppliedFluidMechanics/fluctuations.dat";
    char fluctuationsU_filepath[] = "C:/Users/dlmpa/Numerical_data/5th_semester/AppliedFluidMechanics/fluctuationsU.dat";
    FILE *fluctuations_file,*fluctuationsU_file;
    fluctuations_file = fopen(fluctuations_filepath,"wb");
    fluctuationsU_file = fopen(fluctuationsU_filepath,"wb");

    //percentages of nominal flow rates
    double fluctuations[] = {70 , 75 , 80 , 100 , 120 , 130 , 115 , 105 , 115 , 120 , 120 , 70 };
    scalar_times_vector(12,0.01,fluctuations);

    for(int hour = 0 ; hour < 12 ; hour++){
        //resetting pressures
        double H_reset[] = {2500,
                      2490.0,
                      2485.0,
                      2480.0,
                      2475.0,
                      2470.0,
                      2465.0,
                      2460.0,
                      2455.0,
                      2450.0,
                      2445.0,
                      2440.0,
                      2435.0,
                      2430.0,
                      2425.0,
                      2420.0,
                      2415.0,
                      2410.0,
                      2405.0,
                      2400.0};
        scalar_times_vector(n_NODES,1/0.79/9.81,H_reset);
        //resetting F
        double F_reset[] = {0,0,0,-1.5*0.234 , -1.5*0.234 , -2*0.325 , 0 , -3 * 0.325 , -3 * 0.234 , -2 * 0.234 , -15 * 0.8,
                      0,-1.5*0.234 , -1.5*0.234 , -2*0.325 , 0 , -3 * 0.325 , -3 * 0.234 , -2 * 0.234 , -15 * 0.8};
        scalar_times_vector(n_NODES,fluctuations[hour]*I/3600.0,F_reset);

        newton_solve(n_NODES,n_ELE,ELE,0,F_reset,L,D,Z,K,U,Q,H_reset,1000,1e-4);

        //calculating static pressure
        scalar_times_vector(n_NODES,1*0.79*9.81,H_reset);
        double H_static[n_NODES];
        calc_static(n_NODES,n_ELE,ELE,H_reset,h,U,H_static);

        //saving data to files
        for(int i = 0 ; i < n_NODES ; i++){
            fprintf(fluctuations_file,"%lf ",H_static[i]);
        }
        fprintf(fluctuations_file,"\n");
        for(int i = 0 ; i < n_ELE ; i++){
            fprintf(fluctuationsU_file,"%lf ",U[i]);
        }
        fprintf(fluctuationsU_file,"\n");
    }
    fclose(fluctuations_file);
    fclose(fluctuationsU_file);
    timestop = clock();
    printf("%lf [s] program runtime \n", (double) (timestop-time_start)/CLOCKS_PER_SEC);
    return 0;
}